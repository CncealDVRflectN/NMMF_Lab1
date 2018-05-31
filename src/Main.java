import java.util.Locale;

public class Main {
    private static void printVector(double[] vect) {
        for (double element : vect) {
            System.out.println(element);
        }
    }

    private static void printVectorMath(double[] vect) {
        for (double element : vect) {
            System.out.printf(Locale.ENGLISH, "%E\n", element);
        }
    }

    private static double[] calcDiscrepancy(double[] left, double[] right) {
        double[] result = new double[left.length];

        for (int i = 0; i < result.length; i++) {
            result[i] = left[i] - right[i];
        }

        return result;
    }

    private static double[] calcAlpha(double[] lowerDiagonal, double[] mainDiagonal, double[] upperDiagonal) {
        double[] alpha = new double[mainDiagonal.length - 1];

        alpha[0] = -upperDiagonal[0] / mainDiagonal[0];
        for (int i = 1; i < alpha.length; i++) {
            alpha[i] = -upperDiagonal[i] / (mainDiagonal[i] + lowerDiagonal[i - 1] * alpha[i - 1]);
        }

        return alpha;
    }

    private static double[] calcBeta(double[] lowerDiagonal, double[] mainDiagonal, double[] rightPart, double[] alpha) {
        double[] beta = new double[mainDiagonal.length];

        beta[0] = rightPart[0] / mainDiagonal[0];
        for (int i = 1; i < beta.length; i++) {
            beta[i] = (rightPart[i] - lowerDiagonal[i - 1] * beta[i - 1]) / (mainDiagonal[i] + lowerDiagonal[i - 1] * alpha[i - 1]);
        }

        return beta;
    }

    private static double[] calcRightSweep(double[] lowerDiagonal, double[] mainDiagonal, double[] upperDiagonal, double[] rightPart) {
        double[] result = new double[mainDiagonal.length];
        double[] alpha = calcAlpha(lowerDiagonal, mainDiagonal, upperDiagonal);
        double[] beta = calcBeta(lowerDiagonal, mainDiagonal, rightPart, alpha);

        result[result.length - 1] = beta[beta.length - 1];
        for (int i = result.length - 2; i >= 0; i--) {
            result[i] = alpha[i] * result[i + 1] + beta[i];
        }

        return result;
    }

    /**/

    private static double calcF(double x) {
        return 0.5 - x * x;
    }

    private static double calcP(double x) {
        return -6.0 * x / (3.0 * x * x - 0.5);
    }

    private static double calcQ(double x) {
        return -1.0 / x;
    }

    private static double[] calcFirstOrder(double[] nodes, double step, double[] alpha, double[] beta, double[] gamma) {
        double[] coefsMainDiagonal = new double[nodes.length];
        double[] coefsUpperDiagonal = new double[nodes.length - 1];
        double[] coefsLowerDiagonal = new double[nodes.length - 1];
        double[] rightPart = new double[nodes.length];
        double coef = 1.0 / Math.pow(step, 2.0);
        double pFirst;

        coefsMainDiagonal[0] = alpha[0] - beta[0] / step;
        coefsUpperDiagonal[0] = beta[0] / step;
        rightPart[0] = gamma[0];

        coefsLowerDiagonal[nodes.length - 2] = -beta[1] / step;
        coefsMainDiagonal[nodes.length - 1] = alpha[1] + beta[1] / step;
        rightPart[nodes.length - 1] = gamma[1];

        for (int i = 1; i < nodes.length - 1; i++) {
            pFirst = calcP(nodes[i]) / step;

            coefsLowerDiagonal[i - 1] = coef;
            coefsMainDiagonal[i] = calcQ(nodes[i]) - pFirst - 2.0 * coef;
            coefsUpperDiagonal[i] = coef + pFirst;
            rightPart[i] = calcF(nodes[i]);
        }

        return  calcRightSweep(coefsLowerDiagonal, coefsMainDiagonal, coefsUpperDiagonal, rightPart);
    }

    private static double[] calcAlphaSecondOrder(double[] alpha, double step, double[] betaSecond, double[] nodes) {
        double[] result = new double[alpha.length];

        result[0] = alpha[0] + calcQ(nodes[0]) * betaSecond[0] * step / 2.0;
        result[1] = alpha[1] - calcQ(nodes[nodes.length - 1]) * betaSecond[1] * step / 2.0;

        return result;
    }

    private static double[] calcBetaSecondOrder(double[] beta, double step, double[] nodes) {
        double[] result = new double[beta.length];

        result[0] = beta[0] / (1.0 - calcP(nodes[0]) * step / 2.0);
        result[1] = beta[1] / (1.0 + calcP(nodes[nodes.length - 1]) * step / 2.0);

        return result;
    }

    private static double[] calcGammaSecondOrder(double[] gamma, double step, double[] betaSecond, double[] nodes) {
        double[] result = new double[gamma.length];

        result[0] = gamma[0] + calcF(nodes[0]) * betaSecond[0] * step / 2.0;
        result[1] = gamma[1] - calcF(nodes[nodes.length - 1]) * betaSecond[1] * step / 2.0;

        return result;
    }

    public static double[] calcSecondOrder(double[] nodes, double step, double[] alpha, double[] beta, double[] gamma) {
        double[] coefsMainDiagonal = new double[nodes.length];
        double[] coefsUpperDiagonal = new double[nodes.length - 1];
        double[] coefsLowerDiagonal = new double[nodes.length - 1];
        double[] rightPart = new double[nodes.length];
        double[] betaSecond = calcBetaSecondOrder(beta, step, nodes);
        double[] alphaSecond = calcAlphaSecondOrder(alpha, step, betaSecond, nodes);
        double[] gammaSecond = calcGammaSecondOrder(gamma, step, betaSecond, nodes);
        double coef = 1.0 / Math.pow(step, 2.0);
        double pSecond;

        coefsMainDiagonal[0] = alphaSecond[0] - betaSecond[0] / step;
        coefsUpperDiagonal[0] = betaSecond[0] / step;
        rightPart[0] = gammaSecond[0];

        coefsLowerDiagonal[nodes.length - 2] = - betaSecond[1] / step;
        coefsMainDiagonal[nodes.length - 1] = alphaSecond[1] + betaSecond[1] / step;
        rightPart[nodes.length - 1] = gammaSecond[1];

        for (int i = 1; i < nodes.length - 1; i++) {
            pSecond = calcP(nodes[i]) / (2.0 * step);

            coefsLowerDiagonal[i - 1] = coef - pSecond;
            coefsMainDiagonal[i] = calcQ(nodes[i]) - 2.0 * coef;
            coefsUpperDiagonal[i] = coef + pSecond;
            rightPart[i] = calcF(nodes[i]);
        }

        return calcRightSweep(coefsLowerDiagonal, coefsMainDiagonal, coefsUpperDiagonal, rightPart);
    }

    public static void main(String... args) {
        int splitsNum = 10;
        double intervalBottom = 0.5;
        double intervalUpper = 1.0;
        double step = (intervalUpper - intervalBottom) / splitsNum;
        double[] nodes = new double[splitsNum + 1];
        double[] alpha = new double[2];
        double[] beta = new double[2];
        double[] gamma = new double[2];
        double[] reductionResult = {
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
                0.5053798055106364
        };
        double[] result;

        for (int i = 0; i < nodes.length; i++) {
            nodes[i] = intervalBottom + i * step;
        }

        alpha[0] = 0.0;
        alpha[1] = 2.0;
        beta[0] = 1.0;
        beta[1] = 1.0;
        gamma[0] = 0.25;
        gamma[1] = 3.5;

        System.out.println("Разностная аппроксимация 1-ого порядка:");
        System.out.println("Узлы:");
        printVector(nodes);
        System.out.println("Значения искомой функции в узлax:");
        result = calcFirstOrder(nodes, step, alpha, beta, gamma);
        printVector(result);
        System.out.println("Невязка на решении методом редукции:");
        printVectorMath(calcDiscrepancy(result, reductionResult));
        System.out.println("\nРазностная аппроксимация 2-ого порядка:");
        System.out.println("Узлы:");
        printVector(nodes);
        System.out.println("Значения искомой функции в узлax:");
        result = calcSecondOrder(nodes, step, alpha, beta, gamma);
        printVector(result);
        System.out.println("Невязка на решении методом редукции:");
        printVectorMath(calcDiscrepancy(result, reductionResult));
    }
}
