# Protocol Pas A Pas per a l'Estudi de la Dinàmica Molecular de la LEH

Aquest document proposa un protocol pràctic per estudiar la dinàmica molecular de la LEH amb OpenMM.

Assumpció de treball: `LEH` és l'enzim objecte d'estudi i disposem d'una estructura inicial experimental o modelada.

## 1. Objectiu de l'estudi

Abans de començar, cal definir què volem observar:

- estabilitat global de la LEH
- flexibilitat de llaços o regions actives
- reorganització del centre actiu
- efecte d'un lligand, substrat o mutació
- diferències entre estats apo i holo

Sense aquesta pregunta inicial és fàcil generar trajectòries llargues però poc informatives.

## 2. Preparació de l'estructura inicial

### 2.1. Triar l'estructura de partida

Opcions habituals:

- estructura cristal·logràfica o de cryo-EM
- model d'homologia
- estructura amb substrat o inhibidor
- diferents mutants de la LEH

Criteris recomanats:

- millor resolució possible
- pocs buits estructurals
- presència d'un estat funcional rellevant
- cofactors i ions ben definits

### 2.2. Revisió manual de l'estructura

Cal revisar:

- residus absents
- àtoms faltants
- ocupacions alternatives
- cadenes sobrants
- aigua cristal·logràfica rellevant o prescindible
- lligands, cofactors i ions

### 2.3. Definir l'estat químic

Decisions importants:

- protonació de residus ionitzables a un pH objectiu
- estat de protonació del centre actiu
- estat del substrat o del lligand
- ponts disulfur
- càrrega total del sistema

Per a la LEH això és especialment important si hi ha residus catalítics àcid-base al centre actiu.

## 3. Construcció del sistema

### 3.1. Camp de forces

Recomanació típica:

- proteïna: `amber14-all.xml`
- aigua: `amber14/tip3p.xml`

Si hi ha lligands o substrats:

- cal parametritzar-los amb eines compatibles amb AMBER/OpenMM
- validar càrregues, tipus atòmics i geometria

### 3.2. Afegir hidrògens

Afegir hidrògens segons el pH definit:

- pH fisiològic si no hi ha una justificació diferent
- revisar manualment histidines, glutamats, aspartats, lisines i terminals

### 3.3. Solvatació

Configuració recomanada:

- caixa periòdica
- mínim `1.0-1.2 nm` de padding al voltant de la proteïna
- aigua `TIP3P`

### 3.4. Neutralització i força iònica

Afegir:

- contraions per neutralitzar
- opcionalment `0.15 M` de NaCl per aproximar condicions fisiològiques

## 4. Minimització d'energia

Objectiu:

- eliminar contactes estèrics
- relaxar solvent i ions
- estabilitzar la geometria inicial

Protocol recomanat:

1. minimització amb restriccions suaus sobre la proteïna
2. minimització sense restriccions o amb restriccions molt febles

Valors típics:

- 1000 a 5000 iteracions, o fins convergència

## 5. Equilibratge

L'equilibratge s'ha de fer per etapes.

### 5.1. Etapa NVT

Objectiu:

- escalfar el sistema de manera controlada

Protocol recomanat:

- ensemble NVT
- passar gradualment de `0 K` o `100 K` fins a `300 K`
- restriccions sobre els àtoms pesants de la proteïna
- `50-500 ps`, segons mida i estabilitat

### 5.2. Etapa NPT

Objectiu:

- estabilitzar densitat i volum

Protocol recomanat:

- ensemble NPT a `1 bar`
- mantenir inicialment restriccions sobre la proteïna
- relaxar-les de forma progressiva
- `0.5-2 ns`

### 5.3. Verificacions durant l'equilibratge

Cal revisar:

- temperatura estable
- pressió raonable
- densitat estabilitzada
- absència de col·lapses estructurals
- energia potencial sense derives estranyes

## 6. Simulació de producció

Aquí és on s'obtenen les trajectòries per a l'anàlisi científica.

### 6.1. Condicions recomanades

- ensemble NPT o NVT, segons objectiu
- `2 fs` de timestep amb constriccions sobre enllaços amb H
- PME per electrostàtica de llarg abast
- cutoff no enllaçant de `1.0 nm`

### 6.2. Durada

Per a una LEH:

- exploració inicial: `50-100 ns`
- estudi seriós de flexibilitat: `100-500 ns`
- canvis conformacionals lents: múltiples rèpliques de `200-500 ns` o més

### 6.3. Rèpliques

Molt recomanable fer diverses rèpliques independents:

- mínim 3
- llavors diferents per a velocitats inicials
- mateixes condicions de simulació

Això ajuda a separar soroll estadístic de comportament real.

## 7. Casos d'estudi recomanats per a la LEH

Segons la pregunta biològica, es poden definir diferents campanyes:

### 7.1. LEH apo

Objectiu:

- caracteritzar la flexibilitat basal de l'enzim

### 7.2. LEH amb substrat o lligand

Objectiu:

- veure estabilització del centre actiu
- identificar interaccions persistents
- estudiar l'encaix del substrat

### 7.3. LEH mutant vs wild type

Objectiu:

- comparar estabilitat
- comparar xarxa d'hidrogen
- comparar obertura/tancament de canals

### 7.4. Diferents estats de protonació

Objectiu:

- avaluar sensibilitat del mecanisme a la química del centre actiu

## 8. Anàlisi de trajectòries

L'anàlisi ha d'estar alineada amb la pregunta inicial.

### 8.1. Anàlisi bàsica

- RMSD de la proteïna
- RMSF per residu
- radi de gir
- energia potencial
- nombre de ponts d'hidrogen
- estabilitat de l'estructura secundària

### 8.2. Anàlisi específica del centre actiu

- distàncies entre residus catalítics
- distància substrat-residus clau
- angles rellevants per a la reacció
- presència d'aigua catalítica
- ocupació d'interaccions específiques

### 8.3. Anàlisi de moviments col·lectius

- PCA
- clustering conformacional
- mapes de contactes
- anàlisi de canals o cavitats

### 8.4. Comparacions entre simulacions

Comparar:

- apo vs holo
- wild type vs mutants
- diferents protonacions
- diferents rèpliques

## 9. Criteris de qualitat

Una simulació útil no és només una simulació llarga.

Cal comprovar:

- estabilitat del sistema
- coherència entre rèpliques
- absència d'artefactes evidents
- justificació del model inicial
- consistència entre conclusions i mètriques

## 10. Resultats mínims recomanats per reportar

Per documentar correctament l'estudi:

- estructura inicial utilitzada
- camp de forces
- model d'aigua
- dimensions de la caixa
- nombre total d'àtoms
- temperatura i pressió
- longitud de cada simulació
- nombre de rèpliques
- mètriques analitzades
- figures de RMSD, RMSF i distàncies clau

## 11. Flux de treball resumit

1. obtenir i netejar l'estructura de la LEH
2. definir protonació i estat químic del sistema
3. preparar proteïna, lligands, ions i solvent
4. minimitzar energia
5. equilibrar en NVT
6. equilibrar en NPT
7. executar producció amb diverses rèpliques
8. analitzar estabilitat global
9. analitzar centre actiu i mètriques funcionals
10. comparar condicions o mutants

## 12. Recomanació pràctica per començar

Si vols un protocol senzill però científicament útil, jo començaria així:

1. LEH apo i LEH amb substrat
2. 3 rèpliques per sistema
3. 100 ns per rèplica
4. anàlisi de RMSD, RMSF, distàncies del centre actiu i ponts d'hidrogen

Això et dona una base bona abans de passar a estudis més cars com free energy, QM/MM o EVB.
