# Assemble index.html à partir de index.src.html, solver.js et gpu.js
s = open('index.src.html', encoding='utf-8').read()
s = s.replace(".btn { border:", ".sr { position:absolute; width:1px; height:1px; overflow:hidden; clip:rect(0 0 0 0); }\n.btn { border:", 1)
s = s.replace('/*SOLVER*/', open('solver.js', encoding='utf-8').read())
s = s.replace('/*GPU*/', open('gpu.js', encoding='utf-8').read())
open('index.html', 'w', encoding='utf-8').write(s)
print('index.html généré')
