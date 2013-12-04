import junit.framework.Assert;
import org.junit.Test;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.security.NoSuchAlgorithmException;

public class DiffSolverTest {

    public static final String TEST_PAR_TXT = "testPar.txt";
    public static final String TEST_TXT     = "test.txt";

    public static boolean fileContentsEquals(File file1, File file2) throws IOException {
        FileReader in = null;
        FileReader in1 = null;
        try {
            in = new FileReader(file1);
            in1 = new FileReader(file2);

            BufferedReader br = new BufferedReader(in);
            BufferedReader br2 = new BufferedReader(in1);

            String sCurrentLine;
            String sCurrentLine2;
            int lineNumber = 0;
            while ((sCurrentLine = br.readLine()) != null && (sCurrentLine2 = br2.readLine()) != null) {
                if (!sCurrentLine.equals(sCurrentLine2)) {
                    System.out.println("Line number: " + lineNumber);
                    System.out.println("Simple: " + sCurrentLine);
                    System.out.println("Concurrent: " + sCurrentLine2);
                    return false;
                }
                lineNumber++;
            }
            return true;
        } catch (IOException e) {
            e.printStackTrace();
            throw e;
        } finally {
            try {
                if (in != null) {
                    in.close();
                }
                if (in1 != null) {
                    in1.close();
                }
            } catch (IOException e) {
            }
        }
    }

    @Test
    public void testAllGood() throws IOException, NoSuchAlgorithmException, InterruptedException {
        double u0 = 0.0;
        int f = 1000;
        double dt = 0.001;
        double eps = 0.001;
        int n = 5000;
        int workers = 6;
        DiffSolverParallel solverPar = new DiffSolverParallel(u0, f, dt, eps, n, workers);
        DiffSolver solver = new DiffSolver(u0, f, dt, eps, n);

        solverPar.evolution();
        solver.evolution();

        solverPar.toFile(TEST_PAR_TXT);
        solver.toFile(TEST_TXT);

        Assert.assertTrue("Files are different! ", fileContentsEquals(new File(TEST_TXT), new File(TEST_PAR_TXT)));
    }

    @Test
    public void sequentPerformanceTest() throws IOException {
        double u0 = 0.0;
        int f = 1000;
        double dt = 0.001;
        double eps = 0.001;
        int n = 5000;
        DiffSolver solver = new DiffSolver(u0, f, dt, eps, n);
        long startTime = System.currentTimeMillis();
        solver.evolution();
        System.out.println("Sequent algorithm works: " + (System.currentTimeMillis() - startTime));
    }

    @Test
    public void parallelPerformanceTest() throws IOException, InterruptedException {
        double u0 = 0.0;
        int f = 1000;
        double dt = 0.001;
        double eps = 0.001;
        int n = 5000;
        int workers = 6;
        DiffSolverParallel solverPar = new DiffSolverParallel(u0, f, dt, eps, n, workers);
        long startTime = System.currentTimeMillis();
        solverPar.evolution();
        System.out.println("Parallel algorithm works: " + (System.currentTimeMillis() - startTime));
    }
}
