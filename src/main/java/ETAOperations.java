import Jama.Matrix;

public class ETAOperations {
    public static Matrix get_connectivity(Matrix A, int n, Matrix K1, int reps, int iter) {

        Matrix K2 = A.times(K1);
        Matrix K3;
        double[] dot = MatrixOperations.dotproduct(K2, K1);

        Matrix m1 = new Matrix(iter, reps);
        for (int i = 0; i < K2.getColumnDimension(); i++)
            m1.set(0, i, dot[i]);

        Matrix m2 = new Matrix(iter, reps);
        for (int i = 0; i < iter - 1; i++) {
            Matrix at = MatrixOperations.arraytimes(K1, m1.getMatrix(new int[]{i}, 0, K2.getColumnDimension() - 1));
            K2 = K2.minus(at);
            dot = MatrixOperations.dotproduct(K2, K2);
            Matrix reciprocal = new Matrix(1, reps);
            for (int j = 0; j < K2.getColumnDimension(); j++)
                m2.set(i + 1, j, Math.sqrt(Math.abs(dot[j])));
            for (int j = 0; j < K2.getColumnDimension(); j++)
                reciprocal.set(0, j, 1 / dot[j]);
            K2 = MatrixOperations.arraytimes(K2, reciprocal);
            K3 = A.times(K2).minus(MatrixOperations.arraytimes(K1, m2.getMatrix(new int[]{i + 1}, 0, K2.getColumnDimension() - 1)));
            dot = MatrixOperations.dotproduct(K3, K2);
            for (int j = 0; j < K2.getColumnDimension(); j++)
                m1.set(i + 1, j, dot[i]);
            K1 = K2;
            K2 = K3;
        }

        Matrix ests = new Matrix(1, reps);
        for (int z = 0; z < reps; z++) {
            Matrix x = new Matrix(m2.getRowDimension(), 1);
            for (int j = 0; j < iter - 1; j++) {
                x.set(j, 0, m2.get(j + 1, z));
            }

            Matrix T = MatrixOperations.sparse_matrix(x,
                    m1.getMatrix(0, m1.getRowDimension() - 1, new int[]{z}),
                    m2.getMatrix(0, m2.getRowDimension() - 1, new int[]{z}),
                    iter, iter);
            Matrix U = T.eig().getV();//Get vertices
            Matrix S = T.eig().getD();//Get Dimensions
            //Getting the sparse matrix
            Matrix utemp = U.getMatrix(new int[]{0}, 0, U.getColumnDimension() - 1);
            Matrix m = MatrixOperations.get_diag(MatrixOperations.get_exp(MatrixOperations.get_diag(S)));
            ests.set(0,z, utemp.times(m).times(utemp.transpose()).get(0,0));
        }
        return MatrixOperations.get_log_mean(ests);
    }

    public static double[] equity(Matrix A, Matrix path, int n, double weight, double connectivity_ub, double distance_ub, Matrix b, int reps, int iter, int current_fre, double new_fre, int current_con, double new_con, int full_computation, Matrix base) {
        Matrix A_temp = A;
        double frequency = 0;
        double new_equity = 0;
        double conn = 0;
        Matrix conn_m;
        if(path.getColumnDimension()>1) {
            for (int i = 0; i < path.getColumnDimension() - 2; i++) {
                A_temp.set((int) path.get(i, 0), (int) path.get(i + 1, 0), 1);
                A_temp.set((int) path.get(i + 1, 0), (int) path.get(i, 0), 1);
            }
        }

        conn = current_con + new_con/connectivity_ub;
        if(full_computation==1)
        {
            Matrix temp = get_connectivity(A_temp, n, b, reps, iter).minus(base);
            for(int i = 0;i<temp.getColumnDimension();i++)
                temp.set(0,i,temp.get(0,i)/connectivity_ub);
            conn_m = temp;
        }
        frequency = current_fre + new_fre/distance_ub;
        new_equity = (1 - weight)*conn + weight*frequency;
        return new double[]{frequency, conn, new_equity};
    }

    public static double get_angle(double x0, double y0, double x1, double y1, double x2, double y2) {
        Matrix P0 = new Matrix(new double[][]{{x0, y0}});
        Matrix P1 = new Matrix(new double[][]{{x1, y1}});
        Matrix P2 = new Matrix(new double[][]{{x2, y2}});
        Matrix m1 = P2.minus(P0);
        Matrix m2 = P1.minus(P0);
        double n1 = m1.norm1();
        double n2 = m2.norm1();
        for(int i =0; i<2; i++)
        {
            m1.set(0, i, m1.get(0, i)/n1);
            m2.set(0, i, m1.get(0, i)/n2);
        }
        Matrix n = new Matrix(m1.getRowDimension()+m2.getRowDimension(), m1.getColumnDimension());
        for(int i =0; i<m2.getRowDimension(); i++)
        {
            n.set(i,0,m2.get(i,0));
            n.set(i,1,m2.get(i,1));
        }
        for(int i =0; i<m1.getRowDimension(); i++)
        {
            n.set(i+m2.getRowDimension(),0,m1.get(i,0));
            n.set(i+m2.getRowDimension(),1,m1.get(i,1));
        }
        return (Math.atan2(-1 * n.det(), MatrixOperations.dotproduct(m1,m2)[0]) * 180)/(3.14);
    }
}
