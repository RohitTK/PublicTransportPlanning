import Jama.Matrix;

public class MatrixOperations {
    static double[] dotproduct(Matrix A, Matrix B) {
        Matrix out = new Matrix(1, B.getColumnDimension());
        for (int i = 0; i < A.getColumnDimension(); i++) {
            double sum = 0;
            for (int j = 0; j < B.getRowDimension(); j++) {
                sum += A.get(j, i) * B.get(j, i);
            }
            out.set(0, i, sum);
        }
        return out.getArray()[0];
    }

    static Matrix arraytimes(Matrix A, Matrix B) {
        Matrix out = new Matrix(A.getRowDimension(), B.getColumnDimension());
        for (int i = 0; i < A.getRowDimension(); i++) {
            for (int j = 0; j < A.getColumnDimension(); j++)
                out.set(i, j, A.get(i, j) * B.get(0, j));
        }
        return out;
    }

    public static Matrix sparse_matrix(Matrix l, Matrix m, Matrix r, int row, int col) {
        Matrix out = new Matrix(row, col);
        for (int i = 0; i < m.getRowDimension(); i++) {
            out.set(i, i, m.get(i, 0));
            if (i < col - 1) {
                out.set(i + 1, i, l.get(i, 0));
            }
            if (i > 0) {
                out.set(i - 1, i, r.get(i, 0));
            }
        }
        return out;
    }

    public static Matrix get_diag(Matrix A) {
        Matrix out = new Matrix(A.getRowDimension(), 1);
        if (A.getColumnDimension()==1)
        {
            out = new Matrix(A.getRowDimension(), A.getRowDimension());
            for (int i = 0; i < A.getRowDimension(); i++) {
                out.set(i, 0, A.get(i,0));
            }
        }
        else
        {
            for (int i = 0; i < A.getRowDimension(); i++) {
                out.set(i, 0, A.get(i, i));
            }
        }
        return out;
    }

    public static Matrix get_exp(Matrix A) {
        Matrix out = new Matrix(A.getRowDimension(), A.getColumnDimension());
        for (int i = 0; i < A.getRowDimension(); i++) {
            for (int j = 0; j < A.getColumnDimension(); j++) {
                if (Double.isInfinite(Math.exp(A.get(i, j))))
                    continue;
                out.set(i, j, Math.exp(A.get(i, j)));
            }
        }
        return out;
    }

    public static Matrix get_log_mean(Matrix A) {
        Matrix out = new Matrix(1, A.getColumnDimension());
        for (int i = 0; i < A.getRowDimension(); i++) {
            for (int j = 0; j < A.getColumnDimension(); j++) {
                out.set(0, j, out.get(0, j)+A.get(i, j));
            }
        }
        PublicTransportPlanning.printMatrix(out);
        for (int i = 0; i < A.getColumnDimension(); i++) {
            out.set(0, i, Math.log(Math.abs(out.get(0, i)/A.getRowDimension())));
        }
        return out;
    }
}
