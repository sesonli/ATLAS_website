"""Browser-level responsive checks for every page a visitor can reach.

Started as a named-motif check, which is why only those three routes were
covered when the home page, the user guide and the kink-turn example were all
overflowing a phone viewport. Keep new public routes in ROUTES.
"""

import asyncio
import json
import os
import shutil
import socket
import subprocess
import sys
import tempfile
import time
import urllib.request
from pathlib import Path

import pytest


BASE_DIR = Path(__file__).resolve().parent.parent
ROUTES = (
    "/",
    "/user-guide",
    "/custom-motif-search",
    "/examples",
    "/examples/kink-turn",
    "/examples/sarcin-ricin",
    "/examples/gnra",
)
VIEWPORTS = (1440, 390, 320)


def free_port():
    with socket.socket() as connection:
        connection.bind(("127.0.0.1", 0))
        return connection.getsockname()[1]


def chromium_path():
    configured = os.environ.get("CHROME_PATH")
    candidates = [
        configured,
        shutil.which("chromium"),
        shutil.which("chromium-browser"),
        shutil.which("google-chrome"),
        *sorted(
            Path.home().glob(
                ".cache/ms-playwright/chromium-*/chrome-linux*/chrome"
            ),
            reverse=True,
        ),
    ]
    return next(
        (
            str(path)
            for path in candidates
            if path and Path(path).is_file()
        ),
        None,
    )


def wait_for_url(url, timeout=20):
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        try:
            with urllib.request.urlopen(url, timeout=1) as response:
                if response.status == 200:
                    return
        except Exception:
            time.sleep(0.1)
    raise TimeoutError(f"Timed out waiting for {url}")


async def cdp_request(websocket, request_id, method, params=None):
    await websocket.send(
        json.dumps(
            {
                "id": request_id,
                "method": method,
                "params": params or {},
            }
        )
    )
    while True:
        message = json.loads(await websocket.recv())
        if message.get("id") == request_id:
            if "error" in message:
                raise AssertionError(message["error"])
            return message.get("result", {})


async def audit_layouts(websocket_url, base_url):
    websockets = pytest.importorskip("websockets")
    async with websockets.connect(websocket_url, max_size=2**22) as connection:
        request_id = 0
        for width in VIEWPORTS:
            request_id += 1
            await cdp_request(
                connection,
                request_id,
                "Emulation.setDeviceMetricsOverride",
                {
                    "width": width,
                    "height": 900,
                    "deviceScaleFactor": 1,
                    "mobile": width < 600,
                },
            )
            for route in ROUTES:
                request_id += 1
                await cdp_request(
                    connection,
                    request_id,
                    "Page.navigate",
                    {"url": f"{base_url}{route}"},
                )
                for _ in range(100):
                    request_id += 1
                    state = await cdp_request(
                        connection,
                        request_id,
                        "Runtime.evaluate",
                        {
                            "expression": """
                                ({
                                  readyState: document.readyState,
                                  href: location.href,
                                  h1Count:
                                    document.querySelectorAll('h1').length
                                })
                            """,
                            "returnByValue": True,
                        },
                    )
                    page_state = state["result"].get("value", {})
                    if (
                        page_state.get("readyState") == "complete"
                        and page_state.get("href") == f"{base_url}{route}"
                        and page_state.get("h1Count") == 1
                    ):
                        break
                    await asyncio.sleep(0.05)
                else:
                    raise AssertionError(f"{route} did not finish loading")

                request_id += 1
                result = await cdp_request(
                    connection,
                    request_id,
                    "Runtime.evaluate",
                    {
                        "expression": """
                            (() => ({
                              innerWidth: window.innerWidth,
                              documentWidth:
                                document.documentElement.scrollWidth,
                              h1Count:
                                document.querySelectorAll('h1').length,
                              shellWidth:
                                document.querySelector('.shell')
                                  ?.getBoundingClientRect().width || 0
                            }))()
                        """,
                        "returnByValue": True,
                    },
                )
                layout = result["result"]["value"]
                assert layout["innerWidth"] == width
                assert layout["h1Count"] == 1
                assert layout["shellWidth"] <= width
                assert layout["documentWidth"] <= width + 1, (
                    route,
                    width,
                    layout,
                )


def test_public_pages_have_no_viewport_overflow():
    chrome = chromium_path()
    if not chrome:
        pytest.skip("Chromium is not installed")

    app_port = free_port()
    debug_port = free_port()
    environment = {
        **os.environ,
        "ATLAS_SKIP_STARTUP_CLEANUP": "1",
        "PYTHONUNBUFFERED": "1",
    }
    server = subprocess.Popen(
        [
            sys.executable,
            "-c",
            (
                "from app import app; "
                f"app.run(host='127.0.0.1', port={app_port}, "
                "debug=False, use_reloader=False)"
            ),
        ],
        cwd=BASE_DIR,
        env=environment,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    with tempfile.TemporaryDirectory(prefix="atlas-chrome-") as profile:
        browser = subprocess.Popen(
            [
                chrome,
                "--headless=new",
                "--no-sandbox",
                "--disable-gpu",
                f"--remote-debugging-port={debug_port}",
                f"--user-data-dir={profile}",
                "about:blank",
            ],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        try:
            base_url = f"http://127.0.0.1:{app_port}"
            wait_for_url(f"{base_url}/examples")
            wait_for_url(f"http://127.0.0.1:{debug_port}/json/version")
            with urllib.request.urlopen(
                f"http://127.0.0.1:{debug_port}/json/list"
            ) as response:
                targets = json.load(response)
            websocket_url = targets[0]["webSocketDebuggerUrl"]
            asyncio.run(audit_layouts(websocket_url, base_url))
        finally:
            browser.terminate()
            server.terminate()
            browser.wait(timeout=10)
            server.wait(timeout=10)
