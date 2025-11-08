# save as run_gepia2_headless.py
import asyncio
import os
import re
import yaml
from datetime import datetime
from pathlib import Path
from urllib.parse import quote

from playwright.async_api import async_playwright

# ---- config ----
CONFIG_PATH = "config.yaml"
# tags seen in GEPIA2 docs
TAGS = ["similar", "expdif", "boxplot", "stageplot", "survival", "correlation"]
# CSS candidates where GEPIA2 renders charts (robust to minor changes)
PLOT_SELECTORS = [
    "#plot", ".plot", "#container", "#main", ".echarts", "canvas", "svg"
]
VIEWPORT = {"width": 1400, "height": 1000}
WAIT_SEC = 20

def load_config():
    with open(CONFIG_PATH, "r") as f:
        cfg = yaml.safe_load(f)
    base = cfg["tools"][0]["GEPIA2"]  # e.g., http://...detail.php?gene=$genename&tag=$tagname
    genes = cfg["genes"]
    return base, genes

def build_url(base, gene, tag):
    # ensure safe URL even if gene is like ENSG...
    return base.replace("$genename", quote(gene)).replace("$tagname", quote(tag))

async def capture_plot(page, out_png: Path):
    # try a few selectors; take the first that appears and has non-trivial size
    for sel in PLOT_SELECTORS:
        try:
            el = await page.wait_for_selector(sel, timeout=WAIT_SEC * 1000, state="visible")
            box = await el.bounding_box()
            if box and box["width"] > 200 and box["height"] > 150:
                await el.screenshot(path=str(out_png))
                return True
        except Exception:
            continue
    return False

async def main():
    base, genes = load_config()
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    outdir = Path(f"GEPIA2_outputs_{ts}")
    (outdir / "png").mkdir(parents=True, exist_ok=True)
    (outdir / "html").mkdir(parents=True, exist_ok=True)

    async with async_playwright() as p:
        browser = await p.chromium.launch(headless=True)
        context = await browser.new_context(viewport=VIEWPORT)
        page = await context.new_page()

        for gene in genes:
            for tag in TAGS:
                url = build_url(base, gene, tag)
                png_path = outdir / "png" / f"{gene}_{tag}.png"
                html_path = outdir / "html" / f"{gene}_{tag}.html"

                try:
                    # navigate and wait for network to settle; GEPIA2 uses dynamic JS
                    await page.goto(url, wait_until="networkidle", timeout=WAIT_SEC * 1000)
                    # optional: store the HTML for trace/debug
                    html = await page.content()
                    html_path.write_text(html, encoding="utf-8")

                    ok = await capture_plot(page, png_path)
                    status = "OK" if ok else "NO_PLOT_FOUND"
                except Exception as e:
                    status = f"ERROR: {e.__class__.__name__}"

                print(f"{gene:10s} {tag:12s} -> {status}")

        await context.close()
        await browser.close()

if __name__ == "__main__":
    asyncio.run(main())
