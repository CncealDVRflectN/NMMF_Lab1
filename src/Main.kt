fun calcF(x: Double): Double {
    return 0.5 - Math.pow(x, 2.0)
}

fun calcP(x: Double): Double {
    return -6.0 * x / (3.0 * Math.pow(x, 2.0) - 0.5)
}

fun calcQ(x: Double): Double {
    return -1.0 / x
}

fun calcFirstOrder(nodes: Vector, step: Double, alpha: Vector, beta: Vector, gamma: Vector): Vector {
    val mtr = Matrix(nodes.size, nodes.size)
    val vect = Vector(nodes.size)

    mtr[0][0] = alpha[0] - beta[0] / step
    mtr[0][1] = beta[0] / step
    vect[0] = gamma[0]
    mtr[mtr.lineNum - 1][mtr.columnNum - 2] = -beta[1] / step
    mtr[mtr.lineNum - 1][mtr.columnNum - 1] = alpha[1] + beta[1] / step
    vect[vect.size - 1] = gamma[1]
    for (i in 1 until nodes.size - 1) {
        mtr[i][i - 1] = 1.0 / Math.pow(step, 2.0)
        mtr[i][i] = calcQ(nodes[i]) - calcP(nodes[i]) / step - 2.0 / Math.pow(step, 2.0)
        mtr[i][i + 1] = 1.0 / Math.pow(step, 2.0) + calcP(nodes[i]) / step
        vect[i] = calcF(nodes[i])
    }

    return Matrix.calpRightSweep(mtr, vect)
}

private fun calcBetaSecondOrder(beta: Vector, step: Double, nodes: Vector): Vector {
    val result = Vector(beta.size)

    result[0] = beta[0] / (1.0 - calcP(nodes[0]) * step / 2.0)
    result[1] = beta[1] / (1.0 + calcP(nodes[nodes.size - 1]) * step / 2.0)

    return result
}

private fun calcAlphaSecondOrder(alpha: Vector, step: Double, betaSecond: Vector, nodes: Vector): Vector {
    val result = Vector(alpha.size)

    result[0] = alpha[0] + calcQ(nodes[0]) * betaSecond[0] * step / 2
    result[1] = alpha[1] - calcQ(nodes[nodes.size - 1]) * betaSecond[1] * step / 2

    return result
}

private fun calcGammaSecondOrder(gamma: Vector, step: Double, betaSecond: Vector, nodes: Vector): Vector {
    val result = Vector(gamma.size)

    result[0] = gamma[0] + calcF(nodes[0]) * betaSecond[0] * step / 2
    result[1] = gamma[1] - calcF(nodes[nodes.size - 1]) * betaSecond[1] * step / 2

    return result
}

fun calcSecondOrder(nodes: Vector, step: Double, alpha: Vector, beta: Vector, gamma: Vector): Vector {
    val mtr = Matrix(nodes.size, nodes.size)
    val vect = Vector(nodes.size)
    val betaSecond = calcBetaSecondOrder(beta, step, nodes)
    val alphaSecond = calcAlphaSecondOrder(alpha, step, betaSecond, nodes)
    val gammaSecond = calcGammaSecondOrder(gamma, step, betaSecond, nodes)

    mtr[0][0] = alphaSecond[0] - betaSecond[0] / step
    mtr[0][1] = betaSecond[0] / step
    vect[0] = gammaSecond[0]
    mtr[mtr.lineNum - 1][mtr.columnNum - 2] = -betaSecond[1] / step
    mtr[mtr.lineNum - 1][mtr.columnNum - 1] = alphaSecond[1] + betaSecond[1] / step
    vect[vect.size - 1] = gammaSecond[1]
    for (i in 1 until nodes.size - 1) {
        mtr[i][i - 1] = 1.0 / Math.pow(step, 2.0) - calcP(nodes[i]) / (2.0 * step)
        mtr[i][i] = calcQ(nodes[i]) - 2.0 / Math.pow(step, 2.0)
        mtr[i][i + 1] = 1.0 / Math.pow(step, 2.0) + calcP(nodes[i]) / (2.0 * step)
        vect[i] = calcF(nodes[i])
    }

    return Matrix.calpRightSweep(mtr, vect)
}

fun main(args: Array<String>) {
    val n = 10
    val intervalBottom = 0.5
    val intervalUpper = 1.0
    val step = (intervalUpper - intervalBottom) / n
    val nodes = Vector(n + 1)
    val alpha = Vector(2)
    val beta = Vector(2)
    val gamma = Vector(2)
    val reductionResult = Vector(doubleArrayOf(
            -0.1163045194966492,
            -0.10003278079539082,
            -0.07555464597963049,
            -0.042120829676662225,
            0.0010163171062135407,
            0.05460341500772678,
            0.11938639081301465,
            0.1961106633932895,
            0.2855212512247456,
            0.3883628367057963,
            0.5053798055106364).toTypedArray())
    var result: Vector

    for (i in 0 until nodes.size) {
        nodes[i] = intervalBottom + i * step
    }

    alpha[0] = 0.0
    alpha[1] = 2.0
    beta[0] = 1.0
    beta[1] = 1.0
    gamma[0] = 0.25
    gamma[1] = 3.5

    println("Разностная аппроксимация 1-ого порядка:")
    println("Узлы:")
    nodes.print()
    println("Значения искомой функции в узлax:")
    result = calcFirstOrder(nodes, step, alpha, beta, gamma)
    result.print()
    println("Невязка на решении методом редукции:")
    println((result - reductionResult).norm())
    println()
    println("Разностная аппроксимация 2-ого порядка:")
    println("Узлы:")
    nodes.print()
    println("Значения искомой функции в узлax:")
    result = calcSecondOrder(nodes, step, alpha, beta, gamma)
    result.print()
    println("Невязка на решении методом редукции:")
    println((result - reductionResult).norm())
}