import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class DiffSolverParallel {

    public static final String               TEST_PAR_TXT = "testPar.txt";
    private final       ThreadSafeErrorValue error        = new ThreadSafeErrorValue();
    private final int      NUMBER_OF_WORKERS;
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

    public DiffSolverParallel(double u0, int f, double dt, double eps, int n, int workers) {
        this.u0 = u0;
        this.f = f;
        this.dt = dt;
        this.eps = eps;
        this.n = n;
        this.NUMBER_OF_WORKERS = workers;
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

    public static void main(String[] args) throws IOException, InterruptedException {
        DiffSolverParallel solver = new DiffSolverParallel(0.0, 1000, 0.0001, 0.001, 5000, 6);
        long startTime = System.currentTimeMillis();
        solver.evolution();
        System.out.println(System.currentTimeMillis() - startTime);
        solver.toFile(TEST_PAR_TXT);
    }

    public void evolution() throws InterruptedException {
        uCurrent[n / 2 + n * n / 2] = f;
        double t = 0;
        while (t <= 1 && error.getValue() > eps) {
            List<DiffWorker> threadPool = new ArrayList<DiffWorker>(); // just a collection of created threads
            error.setValue(0.0);
            System.arraycopy(uCurrent, 0, uPrevious, 0, n * n);

            int partSize = n / NUMBER_OF_WORKERS;
            //            long start = System.currentTimeMillis();
            for (int part = 0; part < NUMBER_OF_WORKERS; part++) {

                DiffWorker worker;
                if (part == NUMBER_OF_WORKERS - 1) {
                    worker = new DiffWorker(part * partSize + 1, n - 1, uCurrent, uPrevious);
                } else
                    worker = new DiffWorker(part * partSize + 1, (part + 1) * partSize + 1, uCurrent, uPrevious);

                threadPool.add(worker);
                worker.start();
            }

            for (DiffWorker worker : threadPool) {
                worker.join();
            }
            //            System.out.println("Loops worked: " + (System.currentTimeMillis() - start));

            uCurrent[n / 2 + n * n / 2] += f;
            t += dt;
        }
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

    static class ThreadSafeErrorValue {

        private double value = 1.0;

        public synchronized void add(double value) {
            this.value += value;
        }

        public double getValue() {
            return value;
        }

        public void setValue(Double value) {
            this.value = value;
        }
    }

    class DiffWorker extends Thread {

        private final int      from;
        private final int      to;
        private final double[] uCurrent;
        private final double[] uPrevious;

        DiffWorker(int from, int to, double[] uCurrent, double[] uPrevious) {
            this.from = from;
            this.to = to;
            this.uCurrent = uCurrent;
            this.uPrevious = uPrevious;
        }

        @Override
        public void run() {
            //            long startTime = System.currentTimeMillis();
            double localError = 0.0;
            for (int i = from; i < to; i++) {
                for (int j = 1; j < n - 1; j++) {
                    uCurrent[i + n * j] = a * (uPrevious[i + n * (j + 1)] + uPrevious[i + n * (j - 1)]) +
                                          b * (uPrevious[i + 1 + n * j] + uPrevious[i - 1 + n * j]) +
                                          S[j] * uPrevious[i + n * (j - 1)] + C[i] * uPrevious[i - 1 + n * j] +
                                          (c - S[j] - C[i]) * uPrevious[i + n * j];
                    localError += Math.pow(uCurrent[i + n * j] - uPrevious[i + n * j], 2);
                }
            }
            error.add(localError);
            //            System.out.println(String.format("Thread (%s,%s) ended in %s ms", from, to, System.currentTimeMillis() - startTime));
        }
    }
}