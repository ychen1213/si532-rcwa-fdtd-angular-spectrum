"""Archive the BF/Jones discussion, preserving original text and image bytes.

The earlier 20261106 archives are not rewritten. The source JSONL stays local.
"""
from __future__ import annotations

import argparse
import base64
from collections import Counter
import json
from pathlib import Path

from export_conversation_verbatim import NAMES, local_time, sha, verify_markdown


ROOT = Path(__file__).resolve().parents[1]
FIRST_USER = "你知道CB mark这种alignment mark吗"
LAST_USER = "把以上几段对话原文上传至github"
FULL = "20261107_bf_polarization_verbatim.md"
QA = "20261107_bf_polarization_questions_and_answers_verbatim.md"
PROVENANCE = "20261107_bf_polarization_provenance.json"


def text_of(payload):
    blocks = payload.get("content", [])
    unsupported = [b.get("type") for b in blocks if b.get("type") not in
                   ("input_text", "output_text", "text", "input_image")]
    if unsupported:
        raise ValueError("Unsupported message content: " + str(unsupported))
    return "".join(b.get("text", "") for b in blocks if b.get("type") != "input_image")


def select(source, output):
    raw = source.read_bytes().splitlines(keepends=True)
    records = []
    for line_number, line in enumerate(raw, 1):
        try:
            records.append((line_number, json.loads(line)))
        except json.JSONDecodeError:
            continue
    start = end = None
    for number, record in records:
        p = record.get("payload", {})
        if record.get("type") == "response_item" and p.get("type") == "message" and p.get("role") == "user":
            text = text_of(p)
            if start is None and text.startswith(FIRST_USER):
                start = number
            if start is not None and text.rstrip("\r\n") == LAST_USER:
                end = number
                break
    if start is None or end is None or end <= start:
        raise ValueError("The exact BF discussion range was not found")
    messages = []
    excluded_ui = 0
    for number, record in records:
        if not start <= number <= end or record.get("type") != "response_item":
            continue
        p = record.get("payload", {})
        kind = body = None
        images = []
        blocks = []
        question_count = 0
        if p.get("type") == "message":
            role = p.get("role")
            phase = p.get("phase") or p.get("channel")
            if role == "user":
                body = text_of(p)
                if body.startswith(("<environment_context>", "<external_codex_apps_open_page>")):
                    excluded_ui += 1
                    continue
                kind = "user"
            elif role == "assistant" and phase in ("commentary", "final_answer", "final"):
                body = text_of(p)
                kind = "final_answer" if phase == "final" else phase
            if kind is not None:
                for index, block in enumerate(p.get("content", [])):
                    if block.get("type") == "input_image":
                        url = block["image_url"]
                        if not url.startswith("data:image/png;base64,"):
                            raise ValueError("Expected the original embedded PNG")
                        image_bytes = base64.b64decode(url.split(",", 1)[1], validate=True)
                        rel = "assets/20261106_bf_mark_user_original.png"
                        target = output / rel
                        target.parent.mkdir(parents=True, exist_ok=True)
                        if target.exists() and target.read_bytes() != image_bytes:
                            raise ValueError("An existing image has different bytes")
                        target.write_bytes(image_bytes)
                        images.append(dict(content_index=index, path=rel,
                            bytes=len(image_bytes), sha256=sha(image_bytes),
                            source="Original input_image PNG data URL, decoded without modification"))
                    else:
                        text = block.get("text", "")
                        blocks.append(dict(content_index=index, type=block.get("type"),
                            chars=len(text), utf8_bytes=len(text.encode("utf-8")),
                            sha256=sha(text.encode("utf-8"))))
        elif p.get("type") == "function_call" and p.get("name", "").endswith("request_user_input_async"):
            questions = json.loads(p["arguments"])["questions"]
            rendered = []
            for question in questions:
                text = question["title"]
                if question.get("options"):
                    text += "\n\n选项（按原顺序）：\n" + "\n".join("- " + x for x in question["options"])
                rendered.append(text)
            body = "\n\n".join(rendered)
            kind = "interactive_question"
            question_count = len(questions)
        if kind is not None:
            messages.append(dict(id=f"bf-message-{len(messages) + 1:04d}",
                source_line=number, timestamp_utc=record.get("timestamp"), kind=kind,
                body=body, sha256=sha(body.encode("utf-8")),
                body_chars=len(body), body_utf8_bytes=len(body.encode("utf-8")),
                fenced=bool(images) or body.startswith("<send_user_message_question_reply>"),
                original_text_blocks=blocks, images=images, displayed_question_count=question_count))
    if messages[-1]["body"].rstrip("\r\n") != LAST_USER:
        raise AssertionError("The archive request must be the last archived message")
    return messages, dict(first_source_line=start, cutoff_source_line=end,
        source_prefix_sha256=sha(b"".join(raw[:end])),
        excluded_application_context_messages=excluded_ui)


def render(messages, progress):
    chosen = messages if progress else [m for m in messages if m["kind"] != "commentary"]
    title = "BF mark 与 D4C 偏振矩阵扩展：对话原文"
    if not progress:
        title += "（问答阅读版）"
    counts = Counter(m["kind"] for m in chosen)
    out = [f"# {title}\n\n",
        "> 讨论时间：2026-11-06 至 2026-11-07；归档日期：2026-11-07；时区：Asia/Hong_Kong（UTC+08:00）。\n\n",
        "本文件从本聊天完整会话记录直接提取，范围从最初误称 CB mark 的提问，到本次上传原文的请求。保留更名 BF、几何澄清、用户图片、只讨论 RCWA→D4C 的范围修正、偏振矩阵推导及后续确认。原文不摘要、不改写；公式、代码、措辞、空格、换行与链接照原记录保留。编号、角色、时间及导航为归档附加排版。\n\n",
        "用户交互回复和图片消息的原始文字封装用代码块显示。图片以原始 PNG 字节另存，并在该消息下方展示；每个原始文字块、图片及其内容顺序均在校验清单中记录。交互提问保留工具实际展示的原问题及全部选项，仅添加列表排版。\n\n",
        "早期关于整体几何的描述随后由图片澄清，均按当时原文保留，请按顺序阅读。本次只归档讨论，没有新增 RCWA/FDTD 仿真，也不覆盖此前报告或 2026-11-06 的聊天归档。\n\n",
        f"本版共 {len(chosen)} 条：用户 {counts['user']} 条、助手正式回答 {counts['final_answer']} 条、进度消息 {counts['commentary']} 条、交互澄清 {counts['interactive_question']} 组（共 2 问）。\n\n",
        f"[归档导航](README.md) · [逐条原文与图片校验清单]({PROVENANCE})\n\n",
        "## 正式回答导航\n\n"]
    index = 0
    for m in chosen:
        if m["kind"] == "final_answer":
            index += 1
            out.append(f"- [{local_time(m['timestamp_utc'])} · 第 {index} 次正式回答（原文）](#{m['id']})\n")
    out.append("\n---\n\n")
    for m in chosen:
        out.append(f"<a id=\"{m['id']}\"></a>\n\n")
        out.append(f"## {m['id']} · {NAMES[m['kind']]} · {local_time(m['timestamp_utc'])}\n\n")
        out.append(f"<!-- MESSAGE_BODY_BEGIN:{m['id']} -->\n")
        if m["fenced"]:
            out.append("```xml\n")
        out.append(m["body"])
        if m["fenced"]:
            out.append("\n```")
        out.append(f"\n<!-- MESSAGE_BODY_END:{m['id']} -->\n\n")
        for image in m["images"]:
            out.append("图片附件（原始 PNG；此显示链接为归档附加，原文字未改写）：\n\n")
            out.append(f"![用户提供的 BF mark 图，图4B、408/410]({image['path']})\n\n")
        out.append("---\n\n")
    return "".join(out)


def main(source, output):
    output.mkdir(parents=True, exist_ok=True)
    messages, source_info = select(source, output)
    checks = []
    for name, progress in ((FULL, True), (QA, False)):
        path = output / name
        path.write_bytes(render(messages, progress).encode("utf-8"))
        checks.append(verify_markdown(path, messages, include_progress=progress))
    for m in messages:
        for image in m["images"]:
            if sha((output / image["path"]).read_bytes()) != image["sha256"]:
                raise AssertionError("Original image bytes changed")
    counts = dict(Counter(m["kind"] for m in messages))
    result = dict(scope="BF mark discussion only; source session stays local",
        export_date_hong_kong="2026-11-07", timezone="Asia/Hong_Kong", **source_info,
        message_count=len(messages), counts=counts,
        displayed_question_count=sum(m["displayed_question_count"] for m in messages),
        exactness="All original text blocks are copied in original order without rewriting. Interactive-question titles/options retain original strings. Headers, code fences and usable attachment links are added outside the original text body. Embedded image bytes are decoded unchanged.",
        exclusions="Application environment/page-state, internal instructions/reasoning, tool execution output and agent communication",
        messages=[{k: v for k, v in m.items() if k != "body"} for m in messages],
        markdown_checks=checks, original_attachment_bytes_verified=True)
    (output / PROVENANCE).write_bytes((json.dumps(result, ensure_ascii=False, indent=2) + "\n").encode("utf-8"))
    print(json.dumps(dict(**source_info, counts=counts, files=checks,
        images=[im for m in messages for im in m["images"]]), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-jsonl", required=True, type=Path)
    parser.add_argument("--output", type=Path, default=ROOT / "docs/discussions")
    args = parser.parse_args()
    main(args.source_jsonl, args.output)
