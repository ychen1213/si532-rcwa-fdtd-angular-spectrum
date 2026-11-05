"""Export only user-visible messages from one local chat, without rewriting them.

The private source session is read locally and is never copied to the package.
The cutoff is an exact user message, so later export/publishing progress cannot
silently become part of the archived conversation. Message hashes are checked
against the emitted Markdown before the export is considered complete.
"""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timedelta, timezone
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
FIRST_USER = "完成一个最简单的工作 不要复杂化"
LAST_USER = "把以上对话全部原文上传GitHub 都保留md版本 如果刚才运行的整个程序包没有上传的话也上传"
HK = timezone(timedelta(hours=8))
NAMES = {
    "user": "用户",
    "final_answer": "助手 · 正式回答",
    "commentary": "助手 · 进度消息",
    "interactive_question": "助手 · 交互提问",
}


def sha(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def message_text(payload: dict) -> str:
    blocks = payload.get("content", [])
    if any(block.get("type") not in ("input_text", "output_text", "text") for block in blocks):
        raise ValueError("A message contains non-text material; handle it explicitly")
    return "".join(block.get("text", "") for block in blocks)


def select_messages(source: Path):
    raw_lines = source.read_bytes().splitlines(keepends=True)
    parsed = []
    for line_number, raw in enumerate(raw_lines, 1):
        try:
            parsed.append((line_number, json.loads(raw)))
        except json.JSONDecodeError:
            # A live session may end in a partially written line, after cutoff.
            continue
    start = end = None
    for number, record in parsed:
        p = record.get("payload", {})
        if record.get("type") != "response_item" or p.get("type") != "message" or p.get("role") != "user":
            continue
        body = message_text(p)
        if start is None and body.startswith(FIRST_USER):
            start = number
        if body.rstrip("\r\n") == LAST_USER:
            end = number
    if start is None or end is None or end <= start:
        raise ValueError("The exact initial task and archive-request cutoff were not found")
    selected = []
    excluded_ui = 0
    for number, record in parsed:
        if not start <= number <= end or record.get("type") != "response_item":
            continue
        p = record.get("payload", {})
        kind = None
        body = None
        if p.get("type") == "message":
            role = p.get("role")
            phase = p.get("phase") or p.get("channel")
            if role == "user":
                body = message_text(p)
                if body.startswith(("<environment_context>", "<external_codex_apps_open_page>")):
                    excluded_ui += 1
                    continue
                kind = "user"
            elif role == "assistant" and phase in ("commentary", "final_answer", "final"):
                body = message_text(p)
                kind = "final_answer" if phase == "final" else phase
        elif p.get("type") == "function_call" and p.get("name", "").endswith("request_user_input_async"):
            # These titles/options are displayed to the user by the application.
            # Copy the source strings exactly, adding only list formatting.
            args = json.loads(p["arguments"])
            questions = []
            for question in args["questions"]:
                text = question["title"]
                if question.get("options"):
                    text += "\n\n选项（按原顺序）：\n" + "\n".join("- " + option for option in question["options"])
                questions.append(text)
            body = "\n\n".join(questions)
            kind = "interactive_question"
        if kind is not None:
            selected.append(dict(
                id=f"message-{len(selected) + 1:04d}",
                source_line=number,
                timestamp_utc=record.get("timestamp"),
                kind=kind,
                body=body,
                sha256=sha(body.encode("utf-8")),
                body_utf8_bytes=len(body.encode("utf-8")),
                body_chars=len(body),
                fenced=body.startswith("<send_user_message_question_reply>"),
            ))
    if selected[-1]["body"].rstrip("\r\n") != LAST_USER:
        raise AssertionError("Archive cutoff is not the last selected message")
    return selected, dict(
        first_source_line=start,
        cutoff_source_line=end,
        source_prefix_sha256=sha(b"".join(raw_lines[:end])),
        excluded_application_context_messages=excluded_ui,
    )


def local_time(timestamp: str) -> str:
    return datetime.fromisoformat(timestamp.replace("Z", "+00:00")).astimezone(HK).isoformat(timespec="seconds")


def render(messages: list[dict], *, include_progress: bool) -> str:
    records = messages if include_progress else [m for m in messages if m["kind"] != "commentary"]
    title = "本聊天完整原文记录（含进度消息）" if include_progress else "本聊天问答原文记录（不含进度消息）"
    counts = Counter(m["kind"] for m in records)
    out = [
        f"# {title}\n\n",
        "> 从 2026-10-31 的最初 RCWA/FDTD 任务，到 2026-11-06 的本次原文上传请求。时间采用 Asia/Hong_Kong（UTC+08:00）。\n\n",
        "本文件直接从本聊天的本地会话记录提取，不用摘要或回忆补写。用户输入、助手回答中的措辞、公式、代码、链接、空格及原有更正均保留；只在原文外添加编号、时间和角色。交互提问使用工具向用户展示的原问题及选项文字。用户通过交互界面提交的回复，保留其原始文本封装并用代码块显示。\n\n",
        "原文中可能存在被后续讨论修正的说法，请按时间顺序阅读。原始数学定界符、转义符及本地文件链接未改写；两张原文图片的仓库内链接见 [归档导航](README.md)。本次归档不代表重新运行仿真，也不覆盖原报告。\n\n",
        f"记录数量：{len(records)} 条；用户输入 {counts['user']} 条，助手正式回答 {counts['final_answer']} 条，进度消息 {counts['commentary']} 条，交互提问 {counts['interactive_question']} 条。\n\n",
        "逐条 SHA-256 及导出范围见 [原文校验清单](20261106_conversation_provenance.json)。\n\n",
        "## 正式回答导航\n\n",
    ]
    answer_number = 0
    for m in records:
        if m["kind"] == "final_answer":
            answer_number += 1
            out.append(f"- [{local_time(m['timestamp_utc'])} · 第 {answer_number} 次正式回答（原文）](#{m['id']})\n")
    out.append("\n---\n\n")
    for m in records:
        out.append(f"<a id=\"{m['id']}\"></a>\n\n")
        out.append(f"## {m['id']} · {NAMES[m['kind']]} · {local_time(m['timestamp_utc'])}\n\n")
        out.append(f"<!-- MESSAGE_BODY_BEGIN:{m['id']} -->\n")
        if m["fenced"]:
            out.append("```xml\n")
        out.append(m["body"])
        if m["fenced"]:
            out.append("\n```")
        out.append(f"\n<!-- MESSAGE_BODY_END:{m['id']} -->\n\n---\n\n")
    return "".join(out)


def verify_markdown(path: Path, messages: list[dict], *, include_progress: bool):
    text = path.read_bytes().decode("utf-8")
    expected = messages if include_progress else [m for m in messages if m["kind"] != "commentary"]
    previous = -1
    for m in expected:
        begin = f"<!-- MESSAGE_BODY_BEGIN:{m['id']} -->\n"
        end = f"\n<!-- MESSAGE_BODY_END:{m['id']} -->"
        if text.count(begin) != 1 or text.count(end) != 1:
            raise AssertionError("Missing or duplicate message: " + m["id"])
        start = text.index(begin) + len(begin)
        if start <= previous:
            raise AssertionError("Message order changed")
        stop = text.index(end, start)
        actual = text[start:stop]
        if m["fenced"]:
            if not actual.startswith("```xml\n") or not actual.endswith("\n```"):
                raise AssertionError("Reply display fence changed")
            actual = actual[len("```xml\n"):-len("\n```")]
        if actual != m["body"] or sha(actual.encode("utf-8")) != m["sha256"]:
            raise AssertionError("Verbatim message changed: " + m["id"])
        previous = stop
    if text.count("<!-- MESSAGE_BODY_BEGIN:") != len(expected):
        raise AssertionError("Unexpected additional message body")
    return dict(file=path.name, messages=len(expected), exact_message_text_verified=True,
                sha256=sha(path.read_bytes()), bytes=path.stat().st_size)


def main(source: Path, output: Path):
    messages, source_info = select_messages(source)
    output.mkdir(parents=True, exist_ok=True)
    names = [
        ("20261106_all_messages_verbatim.md", True),
        ("20261106_questions_and_answers_verbatim.md", False),
    ]
    checks = []
    for name, progress in names:
        target = output / name
        target.write_bytes(render(messages, include_progress=progress).encode("utf-8"))
        checks.append(verify_markdown(target, messages, include_progress=progress))
    counts = Counter(m["kind"] for m in messages)
    provenance = dict(
        scope="One user-owned chat, from the initial 532nm task through the exact archive request; source session is not published",
        export_date_hong_kong="2026-11-06",
        timezone="Asia/Hong_Kong",
        **source_info,
        message_count=len(messages),
        counts=dict(counts),
        exclusions="Application environment/page-state envelopes, system/developer instructions, internal reasoning, tool execution records and agent communication are not user-visible conversation",
        exactness="Message bodies are copied without rewriting. Interactive-question titles/options are copied from their displayed source strings, with list formatting added. Markdown reply fences and headings are outside message bodies.",
        no_missing_user_visible_text=True,
        messages=[{k: v for k, v in m.items() if k != "body"} for m in messages],
        markdown_checks=checks,
    )
    (output / "20261106_conversation_provenance.json").write_bytes(
        (json.dumps(provenance, ensure_ascii=False, indent=2) + "\n").encode("utf-8"))
    print(json.dumps(dict(**source_info, counts=dict(counts), files=checks), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-jsonl", required=True, type=Path)
    parser.add_argument("--output", type=Path, default=ROOT / "docs/discussions")
    args = parser.parse_args()
    main(args.source_jsonl, args.output)
