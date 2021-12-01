import Jama.Matrix;

public class Graph {
    public double adjMatrix[][];

    public Graph(int vertex_cnt) {
        adjMatrix = new double[vertex_cnt][vertex_cnt];
    }

    public void addEdge(int v1, int v2, int w) {
        adjMatrix[v1][v2] = adjMatrix[v2][v1] = w;
    }

    public Matrix get_neighbors(int n)
    {
        double[] d = new double[adjMatrix.length];
        int cnt=0;
        for(int i=0;i<adjMatrix.length;i++)
        {
            for(int j=0;j<adjMatrix.length;j++)
            {
                if(adjMatrix[i][j]!=0 && j==n)
                {
                    d[cnt]=adjMatrix[i][j];
                    cnt++;
                }
            }
        }
        Matrix out = new Matrix(cnt, 1);
        for(int i=0;i<cnt;i++)
        {
            out.set(i,0, d[i]);
        }
        return out;
    }

}