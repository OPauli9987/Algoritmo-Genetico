import random
import copy
import time
import matplotlib.pyplot as plt
import blosum

random.seed(42)

blosum62 = blosum.BLOSUM(62)

NFE_ORIGINAL = 0
NFE_MEJORADO = 0


def get_sequences():
    seq1 = "MGSSHHHHHHSSGLVPRGSHMASMTGGQQMGRDLYDDDDKDRWGKLVVLGAVTQGQKLVVLGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQV"

    seq2 = "MKTLLVAAAVVAGGQGQAEKLVKQLEQKAKELQKQLEQKAKELQKQLEQKAKELQKQLEQKAKELQKQLEQKAGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQKELQKQLGQKAKEL"

    seq3 = "MAVTQGQKLVVLGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFAVVAGGQGQAEKLVKQLEQKAKELQKQLEQKAKELQKQLEQKAKELQKQLEQKAKELQKQLEQKALCVFAIN"

    return [list(seq1), list(seq2), list(seq3)]


# FUNCIONES GENERALES

def crear_individuo():
    return get_sequences()


def crear_poblacion_inicial(n=10):
    individuo_base = crear_individuo()
    poblacion = [[row[:] for row in individuo_base] for _ in range(n)]
    return poblacion


def igualar_longitud_secuencias(individuo, gap='-'):
    max_len = max(len(fila) for fila in individuo)
    return [fila + [gap]*(max_len - len(fila)) for fila in individuo]


def validar_poblacion_sin_gaps(poblacion, originales):

    for individuo in poblacion:

        for seq, seq_orig in zip(individuo, originales):

            seq_sin_gaps = [a for a in seq if a != '-']
            seq_orig_sin_gaps = [a for a in seq_orig if a != '-']

            if seq_sin_gaps != seq_orig_sin_gaps:
                return False

    return True


def evaluar_individuo_blosum62(individuo, modo="original"):

    global NFE_ORIGINAL, NFE_MEJORADO

    if modo == "original":
        NFE_ORIGINAL += 1
    else:
        NFE_MEJORADO += 1

    score = 0

    n_seqs = len(individuo)

    seq_len = len(individuo[0])

    for col in range(seq_len):

        for i in range(n_seqs):

            for j in range(i + 1, n_seqs):

                a = individuo[i][col]
                b = individuo[j][col]

                if a == '-' or b == '-':
                    score -= 4
                else:
                    score += blosum62[a][b]

    return score


def eliminar_peores(poblacion, scores, porcentaje=0.5):

    idx_ordenados = sorted(
        range(len(scores)),
        key=lambda i: scores[i],
        reverse=True
    )

    n_seleccionados = max(2, int(len(poblacion) * porcentaje))

    ind_seleccionados = [
        poblacion[i]
        for i in idx_ordenados[:n_seleccionados]
    ]

    scores_seleccionados = [
        scores[i]
        for i in idx_ordenados[:n_seleccionados]
    ]

    return ind_seleccionados, scores_seleccionados


# MUTACIÓN

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


# CRUZA ORIGINAL

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


def cruzar_poblacion_doble_punto(poblacion):

    nueva_poblacion = []

    n = len(poblacion)

    indices = list(range(n))

    random.shuffle(indices)

    parejas = [
        (indices[i], indices[i+1])
        for i in range(0, n-1, 2)
    ]

    if n % 2 == 1:
        parejas.append((indices[-1], indices[0]))

    for idx1, idx2 in parejas:

        padre1 = poblacion[idx1]
        padre2 = poblacion[idx2]

        hijo1, hijo2 = cruzar_individuos_doble_punto(
            padre1,
            padre2
        )

        nueva_poblacion.append(copy.deepcopy(padre1))
        nueva_poblacion.append(copy.deepcopy(padre2))
        nueva_poblacion.append(hijo1)
        nueva_poblacion.append(hijo2)

    return nueva_poblacion[:2*n]


# MEJORAS

# 1. Selección por torneo
# 2. Elitismo
# 3. Mutación adaptativa


def torneo_seleccion(poblacion, scores, k=3):

    candidatos = random.sample(
        list(range(len(poblacion))),
        min(k, len(poblacion))
    )

    ganador = max(candidatos, key=lambda i: scores[i])

    return copy.deepcopy(poblacion[ganador])


def cruzar_individuos_mejorado(
    ind1,
    ind2,
    prob_mutacion=0.25
):

    hijo1 = []
    hijo2 = []

    for seq1, seq2 in zip(ind1, ind2):

        aa1 = [a for a in seq1 if a != '-']
        aa2 = [a for a in seq2 if a != '-']

        if len(aa1) < 8:

            hijo1.append(seq1[:])
            hijo2.append(seq2[:])

            continue

        p1 = random.randint(1, len(aa1)-4)
        p2 = random.randint(p1+2, len(aa1)-1)

        nueva1 = aa1[:p1] + aa2[p1:p2] + aa1[p2:]
        nueva2 = aa2[:p1] + aa1[p1:p2] + aa2[p2:]

        def reconstruir(base_seq, nueva_aa):

            resultado = []

            idx = 0

            for a in base_seq:

                if a == '-':
                    resultado.append('-')
                else:
                    resultado.append(nueva_aa[idx])
                    idx += 1

            return resultado

        hijo1.append(reconstruir(seq1, nueva1))
        hijo2.append(reconstruir(seq2, nueva2))

    hijo1 = mutar_individuo(hijo1, 1, prob_mutacion)
    hijo2 = mutar_individuo(hijo2, 1, prob_mutacion)

    return hijo1, hijo2


def seleccionar_mejorados(
    poblacion,
    scores,
    elite=2
):

    idx_ordenados = sorted(
        range(len(scores)),
        key=lambda i: scores[i],
        reverse=True
    )

    nueva_poblacion = [
        copy.deepcopy(poblacion[i])
        for i in idx_ordenados[:elite]
    ]

    while len(nueva_poblacion) < len(poblacion):

        padre1 = torneo_seleccion(poblacion, scores)
        padre2 = torneo_seleccion(poblacion, scores)

        hijo1, hijo2 = cruzar_individuos_mejorado(
            padre1,
            padre2
        )

        nueva_poblacion.extend([hijo1, hijo2])

    return nueva_poblacion[:len(poblacion)]


# EJECUCIÓN ORIGINAL

def ejecutar_original():

    historial = []

    veryBest = None
    fitnessVeryBest = None

    poblacion = crear_poblacion_inicial(10)

    poblacion = mutar_poblacion_v2(
        poblacion,
        num_gaps=1
    )

    poblacion = [
        igualar_longitud_secuencias(ind)
        for ind in poblacion
    ]

    scores = [
        evaluar_individuo_blosum62(ind, "original")
        for ind in poblacion
    ]

    poblacion, scores = eliminar_peores(
        poblacion,
        scores
    )

    for generaciones in range(30):

        poblacion = cruzar_poblacion_doble_punto(
            poblacion
        )

        poblacion = [
            igualar_longitud_secuencias(ind)
            for ind in poblacion
        ]

        scores = [
            evaluar_individuo_blosum62(ind, "original")
            for ind in poblacion
        ]

        poblacion, scores = eliminar_peores(
            poblacion,
            scores
        )

        best, fitness_best = obtener_best(
            scores,
            poblacion
        )

        if veryBest is None or fitness_best > fitnessVeryBest:

            veryBest = best
            fitnessVeryBest = fitness_best

        historial.append(fitnessVeryBest)

    return historial, fitnessVeryBest, poblacion


# EJECUCIÓN MEJORADA

def ejecutar_mejorado():

    historial = []

    veryBest = None
    fitnessVeryBest = None

    poblacion = crear_poblacion_inicial(14)

    poblacion = mutar_poblacion_v2(
        poblacion,
        num_gaps=2
    )

    poblacion = [
        igualar_longitud_secuencias(ind)
        for ind in poblacion
    ]

    scores = [
        evaluar_individuo_blosum62(ind, "mejorado")
        for ind in poblacion
    ]

    for generaciones in range(30):

        poblacion = seleccionar_mejorados(
            poblacion,
            scores,
            elite=2
        )

        poblacion = [
            igualar_longitud_secuencias(ind)
            for ind in poblacion
        ]

        scores = [
            evaluar_individuo_blosum62(ind, "mejorado")
            for ind in poblacion
        ]

        best, fitness_best = obtener_best(
            scores,
            poblacion
        )

        if veryBest is None or fitness_best > fitnessVeryBest:

            veryBest = best
            fitnessVeryBest = fitness_best

        # Mutación adaptativa
        if generaciones > 5 and len(historial) >= 3:

            if historial[-1] == historial[-2] == historial[-3]:

                for i in range(len(poblacion)):

                    poblacion[i] = mutar_individuo(
                        poblacion[i],
                        2,
                        0.35
                    )

        historial.append(fitnessVeryBest)

    return historial, fitnessVeryBest, poblacion


# OBTENER MEJOR

def obtener_best(scores, poblacion):

    idx_mejor = scores.index(max(scores))

    fitness_best = scores[idx_mejor]

    best = copy.deepcopy(poblacion[idx_mejor])

    return best, fitness_best


if __name__ == "__main__":

    originales = get_sequences()

    inicio = time.time()

    historial_original, fitness_original, poblacion_original = ejecutar_original()

    tiempo_original = time.time() - inicio

    inicio2 = time.time()

    historial_mejorado, fitness_mejorado, poblacion_mejorada = ejecutar_mejorado()

    tiempo_mejorado = time.time() - inicio2

    print("RESULTADOS")

    print("Fitness Original:", fitness_original)
    print("Fitness Mejorado:", fitness_mejorado)

    print("\nTiempo Original:", tiempo_original)
    print("Tiempo Mejorado:", tiempo_mejorado)

    print(
        "\nValidacion Original:",
        validar_poblacion_sin_gaps(
            poblacion_original,
            originales
        )
    )

    print(
        "Validacion Mejorada:",
        validar_poblacion_sin_gaps(
            poblacion_mejorada,
            originales
        )
    )

    # GRÁFICA

    plt.figure(figsize=(12,7))

    plt.plot(
        historial_original,
        marker='o',
        linewidth=3,
        label='Algoritmo Original'
    )

    plt.plot(
        historial_mejorado,
        marker='s',
        linewidth=3,
        label='Algoritmo Mejorado'
    )

    plt.xlabel('Generacion')

    plt.ylabel('Fitness')

    plt.title(
        'Comparacion Algoritmo Original vs Mejorado'
    )

    plt.legend()

    plt.grid(True)

    plt.show()