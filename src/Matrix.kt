class Matrix(var lineNum: Int, var columnNum: Int) {
    var matrix: Array<Array<Double>> = Array(lineNum, { Array(columnNum, { 0.0 }) })

    constructor(initMtr: Array<Array<Double>>) : this(initMtr.size, initMtr[0].size) {
        for (i in matrix.indices) {
            for (j in matrix[0].indices) {
                matrix[i][j] = initMtr[i][j]
            }
        }
    }

    constructor(initMtr: Matrix) : this(initMtr.matrix)

    companion object {
        private fun calcAlpha(mtr: Matrix): Vector {
            val alpha = Vector(mtr.lineNum - 1)

            alpha[0] = -mtr[0][1] / mtr[0][0]
            for (i in 1 until alpha.size) {
                alpha[i] = -mtr[i][i + 1] / (mtr[i][i] + mtr[i][i - 1] * alpha[i - 1])
            }

            return alpha
        }

        private fun calcBeta(mtr: Matrix, vect: Vector, alpha: Vector): Vector {
            val beta = Vector(mtr.lineNum)

            beta[0] = vect[0] / mtr[0][0]
            for (i in 1 until beta.size) {
                beta[i] = (vect[i] - mtr[i][i - 1] * beta[i - 1]) / (mtr[i][i] + mtr[i][i - 1] * alpha[i - 1])
            }

            return beta
        }

        private fun calcSolution(alpha: Vector, beta: Vector): Vector {
            val solution = Vector(beta.size)

            solution[solution.size - 1] = beta[beta.size - 1]
            for (i in solution.size - 2 downTo 0) {
                solution[i] = alpha[i] * solution[i + 1] + beta[i]
            }

            return solution
        }

        private fun isCorrect(mtr: Matrix): Boolean {
            if (mtr.lineNum != mtr.columnNum) {
                return false
            }
            if (Math.abs(mtr[0][0]) < Math.abs(mtr[0][1])) {
                return false
            }
            if (Math.abs(mtr[mtr.lineNum - 1][mtr.columnNum - 1]) < Math.abs(mtr[mtr.lineNum - 1][mtr.columnNum - 2])) {
                return false
            }
            for (i in 1 until mtr.lineNum - 1) {
                if (Math.abs(mtr[i][i]) < (Math.abs(mtr[i][i - 1]) + Math.abs(mtr[i][i + 1]))) {
                    return false
                }
            }
            return true
        }

        fun calpRightSweep(mtr: Matrix, vect: Vector): Vector {
            var alpha: Vector
            var beta: Vector

            if (!isCorrect(mtr)) {
                throw Exception("Incorrect matrix")
            }

            alpha = calcAlpha(mtr)
            beta = calcBeta(mtr, vect, alpha)

            return calcSolution(alpha, beta)
        }
    }

    fun print() {
        for (i in 0 until lineNum) {
            for (j in 0 until columnNum) {
                print(matrix[i][j].toString() + " ")
            }
            println()
        }
    }

    operator fun times(b: Vector): Vector {
        var result: Vector
        if (columnNum != b.size) {
            throw Exception("Incorrect size")
        }
        result = Vector(lineNum)
        for (i in 0 until lineNum) {
            result[i] = (0 until columnNum).sumByDouble { matrix[i][it] * b[it] }
        }
        return result
    }

    operator fun get(i: Int): Array<Double> {
        return matrix[i]
    }
}