"""Regression checks for core report-theme text contrast."""

import re
from pathlib import Path


def _rgb(value):
    value = value.lstrip("#")
    if len(value) == 3:
        value = "".join(channel * 2 for channel in value)
    return tuple(int(value[index : index + 2], 16) for index in (0, 2, 4))


def _luminance(value):
    converted = []
    for channel in _rgb(value):
        channel /= 255
        converted.append(channel / 12.92 if channel <= 0.04045 else ((channel + 0.055) / 1.055) ** 2.4)
    return 0.2126 * converted[0] + 0.7152 * converted[1] + 0.0722 * converted[2]


def _contrast(foreground, background):
    lighter, darker = sorted((_luminance(foreground), _luminance(background)), reverse=True)
    return (lighter + 0.05) / (darker + 0.05)


def test_core_theme_text_meets_wcag_aa():
    page = (Path(__file__).parents[1] / "index.html").read_text(encoding="utf-8")
    roots = re.findall(r":root\{([^}]+)\}", page)
    assert roots, "index.html must define a :root color palette"
    pairs = (
        ("text", "bg"),
        ("muted", "bg"),
        ("accent", "bg"),
        ("text", "panel"),
        ("muted", "panel"),
        ("accent", "panel"),
    )
    for index, root in enumerate(roots):
        mode = "light" if index == 0 else "dark"
        palette = dict(re.findall(r"--([\w-]+):(#[0-9a-fA-F]{3,6})(?:;|$)", root))
        for foreground, background in pairs:
            measured = _contrast(palette[foreground], palette[background])
            assert measured >= 4.5, (
                f"{mode} {foreground}/{background} contrast is {measured:.2f}:1; expected >= 4.5:1"
            )
