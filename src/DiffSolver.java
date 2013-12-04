import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.Arrays;

public class DiffSolver {

    public static final String TEST_TXT = "test.txt";
    private final double   u0;
    private final int      f;
    private final double   dt;
    private final double   eps;
    private final int      n;
    private final double   dx;
    private final double   dy;
    private final double   d;
    private final double   k;
    private final double   a;
    private final double   b;
    private final double   c;
    private final double[] uCurrent;
    private final double[] uPrevious;
    private       double[] C;
    private       double[] S;

    public DiffSolver(double u0, int f, double dt, double eps, int n) {
        this.u0 = u0;
        this.f = f;
        this.dt = dt;
        this.eps = eps;
        this.n = n;
        dx = dy = 1.0 / n;
        d = 10 * dx * dx;
        k = 10 * dy * dy;
        a = k * dt / (dx * dx);
        b = k * dt / (dy * dy);
        c = 2 * dt * d + 1 - 2 * a - 2 * b;
        uCurrent = new double[n * n];
        uPrevious = new double[n * n];
        initStaticData();
    }

    public static void main(String[] args) throws IOException {
        DiffSolver solver = new DiffSolver(0.0, 1000, 0.001, 0.001, 5000);
        long startTime = System.currentTimeMillis();
        solver.evolution();
        System.out.println(System.currentTimeMillis() - startTime);
        solver.toFile(TEST_TXT);
    }

    private void initStaticData() {
        S = new double[n];
        C = new double[n];
        for (int i = 0; i < n; i++) {
            S[i] = 10 * dx * dt * Math.sin(i * dx) / dx;
            C[i] = 10 * dy * dt * Math.cos(i * dy) / dy;
            for (int j = 0; j < n; j++)
                uCurrent[i + n * j] = u0;
        }
    }

    public void evolution() {
        uCurrent[n / 2 + n * n / 2] = f;
        double t = 0;
        double error = 1;
        while ((t <= 1) && (error > eps)) {
            error = 0;
            System.arraycopy(uCurrent, 0, uPrevious, 0, n * n);
            //            long start = System.currentTimeMillis();
            for (int i = 1; i < (n - 1); i++)
                for (int j = 1; j < n - 1; j++) {
                    uCurrent[i + n * j] = a * (uPrevious[i + n * (j + 1)] + uPrevious[i + n * (j - 1)]) +
                                          b * (uPrevious[i + 1 + n * j] + uPrevious[i - 1 + n * j]) +
                                          S[j] * uPrevious[i + n * (j - 1)] + C[i] * uPrevious[i - 1 + n * j] +
                                          (c - S[j] - C[i]) * uPrevious[i + n * j];
                    error = error + Math.pow(uCurrent[i + n * j] - uPrevious[i + n * j], 2);
                }
            //            System.out.println("Loops worked: " + (System.currentTimeMillis() - start));
            uCurrent[n / 2 + n * n / 2] += f;
            t += dt;
        }
    }

    public void toFile(String path) throws IOException {
        File f = new File(path);
        if (!f.exists()) {
            f.createNewFile();
        }
        FileOutputStream fos = new FileOutputStream(f);
        for (int i = 0; i < n; i++) {
            fos.write(Arrays.toString(Arrays.copyOfRange(uCurrent, i * n, i * n + n)).replace("[", "").replace("]", "")
                            .replace(",", " ").replace(".", ",").getBytes());
            fos.write("\n".getBytes());
        }
        fos.close();
    }
}
