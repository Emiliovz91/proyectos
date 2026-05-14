import random
import blosum
import copy
import time
import matplotlib.pyplot as plt

blosum62 = blosum.BLOSUM(62)

def get_sequences():
    seq1 = "MGSSHHHHHHSSGLVPRGSHMASMTGGQQMGRDLYDDDDKDRWGKLVVLGAVTQGQKLVVLGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQV"
    seq2 = "MKTLLVAAAVVAGGQGQAEKLVKQLEQKAKELQKQLEQKAKELQKQLEQKAKELQKQLEQKAKELQKQLEQKAGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQKELQKQLGQKAKEL"
    seq3 = "MAVTQGQKLVVLGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFAVVAGGQGQAEKLVKQLEQKAKELQKQLEQKAKELQKQLEQKAKELQKQLEQKAKELQKQLEQKALCVFAIN"
    return [list(seq1), list(seq2), list(seq3)]

def crear_individuo():
    return get_sequences()

def crear_poblacion_inicial(n=10):
    individuo_base = crear_individuo()
    return [ [row[:] for row in individuo_base] for _ in range(n) ]

def mutar_poblacion_v2(poblacion, num_gaps=1):
    poblacion_mutada = []
    for individuo in poblacion:
        nuevo_individuo = []
        for fila in individuo:
            fila_mutada = fila[:]
            posiciones = set()
            for _ in range(num_gaps):
                pos = random.randint(0, len(fila_mutada))
                while pos in posiciones:
                    pos = random.randint(0, len(fila_mutada))
                posiciones.add(pos)
                fila_mutada.insert(pos, '-')
            nuevo_individuo.append(fila_mutada)
        poblacion_mutada.append(nuevo_individuo)
    return poblacion_mutada

def igualar_longitud_secuencias(individuo, gap='-'):
    max_len = max(len(fila) for fila in individuo)
    return [fila + [gap]*(max_len - len(fila)) for fila in individuo]

def evaluar_individuo_blosum62(individuo, NFE):
    NFE[0] += 1
    score = 0
    n_seqs = len(individuo)
    seq_len = len(individuo[0])
    for col in range(seq_len):
        for i in range(n_seqs):
            for j in range(i+1, n_seqs):
                a = individuo[i][col]
                b = individuo[j][col]
                if a == '-' or b == '-':
                    score -= 4
                else:
                    score += blosum62[a][b]
    return score

def eliminar_peores(poblacion, scores, porcentaje=0.5):
    idx_ordenados = sorted(range(len(scores)), key=lambda i: scores[i], reverse=True)
    n_seleccionados = int(len(poblacion) * porcentaje)
    ind_seleccionados = [poblacion[i] for i in idx_ordenados[:n_seleccionados]]
    scores_seleccionados = [scores[i] for i in idx_ordenados[:n_seleccionados]]
    return ind_seleccionados, scores_seleccionados

def cruzar_individuos_doble_punto(ind1, ind2):
    hijo1 = []
    hijo2 = []
    for seq1, seq2 in zip(ind1, ind2):
        aa_indices = [i for i, a in enumerate(seq1) if a != '-']
        if len(aa_indices) < 6:
            hijo1.append(seq1[:])
            hijo2.append(seq2[:])
            continue
        intentos = 0
        while True:
            p1, p2 = sorted(random.sample(aa_indices, 2))
            if p2 - p1 >= 5 or intentos > 10:
                break
            intentos += 1

        def cruza(seqA, seqB):
            aaA = [a for a in seqA if a != '-']
            aaB = [a for a in seqB if a != '-']
            nueva = aaA[:p1] + aaB[p1:p2] + aaA[p2:]
            resultado = []
            idx = 0
            for a in seqA:
                if a == '-':
                    resultado.append('-')
                else:
                    resultado.append(nueva[idx])
                    idx += 1
            return resultado

        nueva_seq1 = cruza(seq1, seq2)
        nueva_seq2 = cruza(seq2, seq1)
        hijo1.append(nueva_seq1)
        hijo2.append(nueva_seq2)
    hijo1 = mutar_individuo(hijo1, 1, 0.8)
    hijo2 = mutar_individuo(hijo2, 1, 0.8)
    return hijo1, hijo2

def mutar_individuo(individuo, n_gaps, p):
    nuevo_individuo = []
    for secuencia in individuo:
        sec = secuencia[:]
        if random.random() < p:
            posiciones = set()
            for _ in range(n_gaps):
                pos = random.randint(0, len(sec))
                while pos in posiciones:
                    pos = random.randint(0, len(sec))
                posiciones.add(pos)
                sec.insert(pos, '-')
        nuevo_individuo.append(sec)
    return nuevo_individuo

def cruzar_poblacion_doble_punto(poblacion):
    nueva_poblacion = []
    n = len(poblacion)
    indices = list(range(n))
    random.shuffle(indices)
    parejas = [(indices[i], indices[i+1]) for i in range(0, n-1, 2)]
    if n % 2 == 1:
        parejas.append((indices[-1], indices[0]))
    for idx1, idx2 in parejas:
        padre1 = poblacion[idx1]
        padre2 = poblacion[idx2]
        hijo1, hijo2 = cruzar_individuos_doble_punto(padre1, padre2)
        nueva_poblacion.append(copy.deepcopy(padre1))
        nueva_poblacion.append(copy.deepcopy(padre2))
        nueva_poblacion.append(hijo1)
        nueva_poblacion.append(hijo2)
    return nueva_poblacion[:2*n]

def validar_poblacion_sin_gaps(poblacion, originales):
    for individuo in poblacion:
        for seq, seq_orig in zip(individuo, originales):
            seq_sin_gaps = [a for a in seq if a != '-']
            seq_orig_sin_gaps = [a for a in seq_orig if a != '-']
            if seq_sin_gaps != seq_orig_sin_gaps:
                return False
    return True

def obtener_best(scores, poblacion):
    idx_mejor = scores.index(max(scores))
    fitness_best = scores[idx_mejor]
    best = copy.deepcopy(poblacion[idx_mejor])
    return best, fitness_best

# NUEVO: refrescar población
def refrescar_poblacion(poblacion, num_nuevos):
    for _ in range(num_nuevos):
        individuo_nuevo = crear_individuo()
        poblacion.append([row[:] for row in individuo_nuevo])
    return poblacion

def ejecutar_esquema(nombre, tamanio_poblacion, gaps_mutacion, generaciones, uso_refresco=False, freq_refresco=5, nuevos_refresco=4):
    NFE = [0]
    veryBest = None
    fitnessVeryBest = None
    poblacion = crear_poblacion_inicial(tamanio_poblacion)
    poblacion = mutar_poblacion_v2(poblacion, num_gaps=gaps_mutacion)
    poblacion = [igualar_longitud_secuencias(ind) for ind in poblacion]
    scores = [evaluar_individuo_blosum62(ind, NFE) for ind in poblacion]
    poblacion, scores = eliminar_peores(poblacion, scores)
    best_fitness_por_gen = []
    for gen in range(generaciones):
        poblacion = cruzar_poblacion_doble_punto(poblacion)
        if uso_refresco and (gen+1)%freq_refresco == 0:
            poblacion = refrescar_poblacion(poblacion, nuevos_refresco)
        poblacion = [igualar_longitud_secuencias(ind) for ind in poblacion]
        scores = [evaluar_individuo_blosum62(ind, NFE) for ind in poblacion]
        poblacion, scores = eliminar_peores(poblacion, scores)
        best, fitness_best = obtener_best(scores, poblacion)
        if veryBest is None or fitness_best>fitnessVeryBest:
            veryBest = best
            fitnessVeryBest = fitness_best
        best_fitness_por_gen.append(fitnessVeryBest)
    print(f"\n{'='*50}")
    print(f"Esquema: {nombre}")
    print("Fitness final:", fitnessVeryBest)
    print ("Validacion integridad:", validar_poblacion_sin_gaps(poblacion,get_sequences()))
    print("="*50)
    return best_fitness_por_gen

if __name__ == "__main__":
    generaciones = 30
    # Código original (sin refresco)
    fit_sin_refresco = ejecutar_esquema(
        "original (sin refresco)",
        tamanio_poblacion=15,
        gaps_mutacion=2,
        generaciones=generaciones,
        uso_refresco=False
    )
    # Código mejorado (con refresco)
    fit_con_refresco = ejecutar_esquema(
        "mejorado (con refresco cada 5 gens)",
        tamanio_poblacion=15,
        gaps_mutacion=2,
        generaciones=generaciones,
        uso_refresco=True,
        freq_refresco=5,
        nuevos_refresco=4
    )
    # Graficar
    plt.figure()
    plt.plot(range(generaciones), fit_sin_refresco, label='Sin refresco', color='red', linestyle='--')
    plt.plot(range(generaciones), fit_con_refresco, label='Con refresco', color='green', linestyle='-')
    plt.xlabel('Generación')
    plt.ylabel('Mejor fitness')
    plt.title('Comparativa: Algoritmo Genético\nsin refresco VS con refresco de población')
    plt.legend()
    plt.show()