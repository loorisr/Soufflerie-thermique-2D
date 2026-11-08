import asyncio, json, sys
from playwright.async_api import async_playwright
async def main():
    async with async_playwright() as p:
        b = await p.chromium.launch(args=["--enable-unsafe-webgpu","--enable-features=Vulkan","--use-vulkan=swiftshader","--use-webgpu-adapter=swiftshader","--ignore-gpu-blocklist"])
        pg = await b.new_page(); errs=[]
        pg.on('pageerror', lambda e: errs.append(str(e))); pg.on('console', lambda m: errs.append(m.text))
        await pg.goto('http://localhost:8765/tests/cmp.html')
        for name, steps in [('pac', 1500), ('piece', 1500)]:
            try:
                r = await pg.evaluate(f"runCmp('{name}', {steps})")
                print(name, json.dumps(r, indent=None))
            except Exception as e: print(name, 'ERR', str(e)[:800])
        print('console', errs[:10]); await b.close()
asyncio.run(main())
