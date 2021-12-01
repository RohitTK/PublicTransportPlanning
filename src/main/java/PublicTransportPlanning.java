import Jama.Matrix;
import tech.tablesaw.api.Table;

import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

public class PublicTransportPlanning {
    //Datasets
    static Matrix neighbor_stops_data;
    static Matrix stop_locations_data;
    static Matrix weighted_data;
    static Matrix k1_data;

    static Matrix objectives;

    //parameters
    static int nodesNumber = 12340;
    static double weight = 0.5;
    static int k = 30;
    static int number_of_turns = 1;
    static int seeding_number = 5000;

    static int max_len = k;
    static int max_iter = 100000;
    static int area_label = 5; // 1: Brooklyn, 2: Manhattan, 3: Staten, 4: Queens, 5 Bronx
    static int area_indicator = 1;
    static int city_choice = 1;

    static Matrix table_to_matrix(Table t) {
        Matrix m = new Matrix(t.rowCount(), t.columnCount());
        for (int i = 0; i < t.rowCount(); i++) {
            for (int j = 0; j < t.columnCount(); j++) {
                m.set(i, j, Double.parseDouble(t.get(i, j).toString()));
            }
        }
        return m;
    }

    //Initializing the datasets
    static void init_datasets(String city) throws IOException {
        if (Objects.equals(city, "Chicago")) {
            //Reading Chicago Data
            area_indicator = 0;
            Table neighbor_stops_data_tab = Table.read().file(System.getProperty("user.dir") + "\\Datasets\\chi_data_csv\\chi_transit_neighbors0.5_0.2.csv");
            neighbor_stops_data = table_to_matrix(neighbor_stops_data_tab);
            Table stop_locations_data_tab = Table.read().file(System.getProperty("user.dir") + "\\Datasets\\chi_data_csv\\chi_transit_stop_locations.csv");
            stop_locations_data = table_to_matrix(stop_locations_data_tab);
            Table weighted_data_tab = Table.read().file(System.getProperty("user.dir") + "\\Datasets\\chi_data_csv\\chi_transit_weighted.csv");
            weighted_data = table_to_matrix(weighted_data_tab);
            Table k1_data_tab = Table.read().file(System.getProperty("user.dir") + "\\Datasets\\chi_data_csv\\new_chi_random_K1.csv");
            k1_data = table_to_matrix(k1_data_tab);
        } else {
            //Reading NYC Data
            area_indicator = 1;
            Table neighbor_stops_data_tab = Table.read().file(System.getProperty("user.dir") + "\\Datasets\\nyc_data_csv\\new_nyc_transit_neighbors0.5_0.2_pre.csv");
            neighbor_stops_data = table_to_matrix(neighbor_stops_data_tab);
            Table stop_locations_data_tab = Table.read().file(System.getProperty("user.dir") + "\\Datasets\\nyc_data_csv\\nyc_transit_stop_locations.csv");
            stop_locations_data = table_to_matrix(stop_locations_data_tab);
            Table weighted_data_tab = Table.read().file(System.getProperty("user.dir") + "\\Datasets\\nyc_data_csv\\nyc_transit_weighted.csv");
            weighted_data = table_to_matrix(weighted_data_tab);
            Table k1_data_tab = Table.read().file(System.getProperty("user.dir") + "\\Datasets\\nyc_data_csv\\new_nyc_random_K1.csv");
            k1_data = table_to_matrix(k1_data_tab);
        }
    }

    public static void printMatrix(Matrix m) {
        for (int r = 0; r < m.getRowDimension(); r++) {
            for (int c = 0; c < m.getColumnDimension(); c++) {
                System.out.print(m.get(r, c));
                System.out.print("\t");
            }
            System.out.println();
        }
    }

    public static void main(String[] args) throws IOException {
        int full_computation = 0;
        String city = "Chicago";
        init_datasets("Chicago");
        Matrix B = neighbor_stops_data;
        Hashtable<Integer, Integer> neighbors_node = new Hashtable<Integer, Integer>();
        Hashtable<Integer, Integer> neighbor_number = new Hashtable<Integer, Integer>();

        int A = 0;
        int cursor = 0;
        double[] ranked_weight;

        int i;
        for (i = 0; i < B.getColumnDimension(); i++) {
            if (Math.abs((int) B.get(i, 0) - A) < 1) { //Approximately enough
                neighbors_node.put((int) B.get(i, 0), i);
                if (i > 0) {
                    neighbor_number.put((int) B.get(i - 1, 0), i - cursor - 1);
                }
                cursor = i;
            }
            A = (int) B.get(i, 0);
        }
        neighbor_number.put((int) B.get(i - 1, 0), i - cursor);

        Matrix C = new Matrix(weighted_data.getRowDimension(), weighted_data.getColumnDimension() + 3);

        for (i = 0; i < weighted_data.getRowDimension(); i++)
            for (int j = 0; j < weighted_data.getColumnDimension(); j++)
                C.set(i, j, weighted_data.get(i, j));

        Graph G = new Graph(C.getRowDimension());

//        List<Integer> arr_values = new ArrayList<Integer>();
        for (i = 0; i < C.getRowDimension(); i++) {
            if ((int) C.get(i, 0) > nodesNumber || (int) C.get(i, 1) > nodesNumber) {
                continue;
            }
            int n1 = (int) C.get(i, 0);
            int n2 = (int) C.get(i, 1);
            G.addEdge(n1, n2, i + 1);
        }
        Matrix Adj = new Matrix(G.adjMatrix);

        int graph_dimension = Adj.getRowDimension();

        Matrix D = new Matrix(B.getRowDimension() + C.getRowDimension(), C.getColumnDimension());
        for (i = 0; i < B.getRowDimension(); i++)
            for (int j = 0; j < B.getColumnDimension(); j++)
                D.set(i, j, B.get(i, j));
        for (i = B.getRowDimension(); i < D.getRowDimension(); i++)
            for (int j = 0; j < D.getColumnDimension(); j++)
                D.set(i, j, C.get(i - B.getRowDimension(), j));

        int reps = 50;
        int iter = 10;

        Matrix K1 = k1_data;
        Matrix base = ETAOperations.get_connectivity(Adj, graph_dimension, K1, reps, iter);

        double[][] temp = D.getMatrix(0, D.getRowDimension() - 1, new int[]{5}).getArray();
        double[][] sorttemp = new double[D.getRowDimension()][2];
        for (int j = 0; j < D.getRowDimension(); j++) {
            sorttemp[j][0] = j;
            sorttemp[j][1] = temp[j][0];
        }
        java.util.Arrays.sort(sorttemp, (b, a) -> Double.compare(a[1], b[1]));
        Matrix idxA = new Matrix(sorttemp).getMatrix(0, k, new int[]{0});
        Matrix SortA = new Matrix(sorttemp).getMatrix(0, k, new int[]{1});

        double lb_weight_max = 0;
        for (int j = 0; j < max_len; j++) {
            lb_weight_max += SortA.get(j, 0);
        }

        temp = B.getMatrix(0, B.getRowDimension() - 1, new int[]{6}).getArray();
        sorttemp = new double[temp.length][2];
        for (int j = 0; j < temp.length; j++) {
            sorttemp[j][0] = j;
            sorttemp[j][1] = temp[j][0];
        }
        java.util.Arrays.sort(sorttemp, (b, a) -> Double.compare(a[1], b[1]));

        Matrix idx = new Matrix(sorttemp).getMatrix(0, k, new int[]{0});
        Matrix SortB = new Matrix(sorttemp).getMatrix(0, k, new int[]{1});
        Matrix bound_conn = SortB.getMatrix(0, k, new int[]{0});
        double bottom_conn = SortB.get(k, 0);
        double connectivity_ub = 0;

        for (int j = 0; j < k; j++) {
            if (((Double) SortB.get(j, 0)).isNaN())
                continue;
            connectivity_ub += SortB.get(j, 0);
        }
        double lb_conn_max = connectivity_ub;

        for (int j = 0; j < D.getRowDimension(); j++) {
            D.set(j, 7,
                    (weight * D.get(j, 5) / lb_weight_max)
                            + (1 - weight) * D.get(j, 6) / lb_conn_max);
        }

        temp = D.getMatrix(0, D.getRowDimension() - 1, new int[]{7}).getArray();
        sorttemp = new double[D.getRowDimension()][2];
        for (int j = 0; j < D.getRowDimension(); j++) {
            sorttemp[j][0] = j;
            sorttemp[j][1] = temp[j][0];
        }
        java.util.Arrays.sort(sorttemp, (b, a) -> Double.compare(a[1], b[1]));
        idxA = new Matrix(sorttemp).getMatrix(0, sorttemp.length - 1, new int[]{0});
        SortA = new Matrix(sorttemp).getMatrix(0, sorttemp.length - 1, new int[]{1});
        Matrix bound_weight = SortA.getMatrix(0, max_len - 1, new int[]{0});
        double bottom_weight = SortB.get(max_len - 1, 0);

        double estimated_running_time = 0;
        double time_cost_conn_compute = 0.09;
        double time_cost_ub_compute = 0.13;

        if (area_indicator == 0) {
            time_cost_conn_compute = 0.05;
            time_cost_ub_compute = 0.08;
        }

        int optimization_choice = 0;
        if (full_computation == 1)
            optimization_choice = 0;

        System.out.println("running experiments on: " + k + "_" + weight + "_" + number_of_turns + "_" + seeding_number + "_" + optimization_choice + "_" + full_computation + "_" + city_choice);

        Matrix best_path = null;// = new Matrix(1, 2);
        double min_equity = 0;
        double n = 1;
        PriorityQueue<double[]> q = new PriorityQueue<double[]>(Comparator.comparingDouble(o -> o[1]));
        double lb = 1;
        Hashtable<Integer, Integer> checked_edges = new Hashtable<Integer, Integer>();
        Hashtable<String, Integer> seeding_edges = new Hashtable<String, Integer>();
        Matrix newpath;
        for (i = 0; i < seeding_number; i++) {
            double ix = idxA.get(i, 0);
            newpath = new Matrix(new double[][]{{D.get((int) ix, 0), D.get((int) ix, 1)}});

            if ((seeding_edges.containsKey(D.get((int) ix, 0) + "," + D.get((int) ix, 1)))
                    || (seeding_edges.containsKey(D.get((int) ix, 1) + "," + D.get((int) ix, 0))))
                continue;

            seeding_edges.put((D.get((int) ix, 0) + "," + D.get((int) ix, 1)), 1);
            seeding_edges.put((D.get((int) ix, 1) + "," + D.get((int) ix, 0)), 1);

            double bw = max_len;
            if (D.get((int) ix, 7) < bottom_weight) {
                bw -= 1;
                lb = 1 - (bottom_weight - D.get((int) ix, 7));
            }

            double a[] = ETAOperations.equity(Adj, newpath, graph_dimension, weight, connectivity_ub, lb_weight_max, K1, reps, iter, 0, D.get((int) ix, 5), 0, D.get((int) ix, 6), 0, base);
            double frequency = a[0];
            double conn = a[1];
            double new_equity = a[2];

            if (new_equity > min_equity) {
                min_equity = new_equity;
                best_path = newpath;
            }

            if (area_indicator == 1) {
                if (Math.abs(area_label - stop_locations_data.get((int) D.get((int) ix, 0), 4)) == 1
                        || Math.abs(area_label - stop_locations_data.get((int) D.get((int) ix, 1), 4)) == 1)// Approximately equal
                    continue;
            }
            q.add(new double[]{lb * -1, bw, 0, frequency, conn, (double) D.get((int) ix, 0), (double) D.get((int) ix, 1)});
        }

        int bad = 0;
        double[] candidate;
        while (n < max_iter && q.size() > 0) {
            candidate = q.remove();
            double lower_bound = candidate[0] * -1;

            if (min_equity > lower_bound) {
                break;
            }
            double size_path = candidate.length; //Path size
            double bw = candidate[1];
            double bc = candidate[2];
            double frequency = candidate[3];
            double conn = candidate[4];
            double[] path = Arrays.copyOfRange(candidate, 5, (int) size_path-1);

            if (n % 100 == 0) {
                double cc = n / 100;
                if (city.equals("NYC"))
                    objectives.set((int) cc, optimization_choice + 1 + 3, min_equity);
            }
            double last_node = candidate[(int) size_path-1];
            double first_node = candidate[5];

            Matrix N = G.get_neighbors((int) last_node).transpose();
            List<Double> ne = new ArrayList<Double>();
            for (i = 0; i < N.getColumnDimension()-1; i++) {
                double d = G.adjMatrix[(int) last_node][(int) N.get(0, i)];
                d += B.getRowDimension();
                ne.add(d);
            }
            if (neighbors_node.containsKey(String.valueOf(first_node))) {
                for (i = 0; i < neighbor_number.get(String.valueOf(first_node)); i++) {
                    ne.add((double) (neighbors_node.get(String.valueOf(first_node)) + i));
                }
            }

            if (neighbors_node.containsKey(String.valueOf(last_node))) {
                for (i = 0; i < neighbor_number.get(String.valueOf(last_node)); i++) {
                    ne.add((double) (neighbors_node.get(String.valueOf(last_node)) + i));
                }
            }
            int max_increment = 0;
            int best_xx = 0;
            double[] array = new double[0];
            double bcc = bc;
            int best_neighbor = 0;
            for (i = 0; i < ne.size(); i++) {
                double index_neighbor = ne.get(i);
                double xx = D.get((int) index_neighbor - 1, 1);
                if (xx == last_node) {
                    xx = D.get((int) index_neighbor - 1, 0);
                }
                if (area_indicator == 1) {
                    if (Math.abs(area_label - stop_locations_data.get((int) xx - 1, 3)) < 1) {
                        continue;
                    }
                }
                for (double d : path) {
                    if (xx == d)
                        continue;
                }

                double new_bc = bc;
                double last_s_node = candidate[(int) (size_path - 2)];
                double angle = Math.abs(ETAOperations.get_angle(stop_locations_data.get((int) last_node, 2),
                        stop_locations_data.get((int) last_node, 1),
                        stop_locations_data.get((int) xx, 2),
                        stop_locations_data.get((int) xx, 1),
                        stop_locations_data.get((int) last_s_node, 2),
                        stop_locations_data.get((int) last_s_node, 1)));

                angle = Double.parseDouble(String.valueOf(angle).substring(0,3));
                if (angle < (3.14)/2)
                    continue;
                if (angle < 3*(3.14)/4)
                    new_bc = bc + 1;
                if (new_bc > number_of_turns)
                    continue;


                newpath = new Matrix(1, path.length+1);
                printMatrix(newpath);
                System.out.println(xx);
                for(i=newpath.getColumnDimension()-1;i>1; i--)
                    newpath.set(0,i, newpath.get(0,i));
                newpath.set(0,0,xx);
                printMatrix(newpath);

                double[] temp_e = ETAOperations.equity(Adj, newpath, graph_dimension, weight, connectivity_ub, lb_weight_max, K1, reps, iter, (int)frequency, D.get((int) index_neighbor, 6), (int) conn, D.get((int) index_neighbor, 7), 0, base);
                double newfre = temp_e[0];
                double newconn = temp_e[1];
                double new_equity = temp_e[1];
                estimated_running_time = estimated_running_time + time_cost_conn_compute;

                double[] bestarray = new double[10];

                if (new_equity > min_equity) {
                    min_equity = new_equity;
                    best_path = newpath;
                    bestarray = candidate;
                    bestarray[4] = newfre;
                    bestarray[5] = newconn;
                }


                double incrementcc = D.get((int) index_neighbor, 7);
                if (incrementcc > max_increment) {
                    max_increment = (int) incrementcc;
                    lb = lower_bound;
                    bottom_weight = bound_weight.get(0, (int) bw);
                    double new_bw=0;
                    if(D.get((int) index_neighbor,7)<bottom_weight)
                    {
                        new_bw = bw-1;
                        lb = lb - (bottom_weight-D.get((int) index_neighbor,7));
                    }
                    if (new_bw >=1)
                    {
                        double[] t = new double[5 + newpath.getColumnDimension()];
                        t[0] = lb * -1;
                        t[1] = new_bw;
                        t[2] = new_bc;
                        t[3] = newfre;
                        t[4] = newconn;
                        for(i=5; i<t.length; i++)
                        {
                            t[i]=newpath.get(0 ,i-5);
                        }
                        array = t;
                    }
                    if (optimization_choice == 1)
                    {
                        if ((array.length < (5 + max_len)) && (array.length >=7))
                            q.add(array);
                        estimated_running_time = estimated_running_time + time_cost_ub_compute;
                    }
                    else
                    {
                        bad += 1;
                        int a=0;
                    }
                    n=n+1;
                }

                int count_new_edges = 0;

                System.out.println("#Bad: "+ bad);
                System.out.println("We found the path: ");
                for(i=0;i<best_path.getColumnDimension();i++){
                    System.out.print(best_path.get(0,i)+" ");
                }
                
                FileWriter fileIDnew = new FileWriter("./new_edge_b.txt");
                Matrix A_temp = Adj;
                String content = "geometry\\n\"{\"\"type\"\": \"\"LineString\"\", \"\"coordinates\"\": [";
                double id;
                for(i=0;i<best_path.getColumnDimension()-1;i++) {
                    if(Adj.get((int) best_path.get(0,i), (int) best_path.get(0, i+1)) == 0)
                    {
                        id = best_path.get(0,i);
                        count_new_edges+=1;
                        fileIDnew.write(String.valueOf(stop_locations_data.get((int)id,2)) + String.valueOf(stop_locations_data.get((int)id,1)));
                    }
                    A_temp.set((int)best_path.get(0, i),(int)best_path.get(0, i+1), 1);
                    A_temp.set((int)best_path.get(0, i+1),(int)best_path.get(0, i), 1);

                }

                System.out.println(bestarray.toString());
                System.out.println("Increased Equity: " + min_equity);
                System.out.println("#New Edges: " + count_new_edges);
                Matrix base1 = ETAOperations.get_connectivity(A_temp,graph_dimension, K1, reps, iter);

                System.out.println("Estimated Increase Conn: "+bestarray[4]);

                System.out.println("Estimated increase demand: "+bestarray[3]);
                System.out.println("Estimated time on conn and ub: "+ estimated_running_time);

                FileWriter fileID = new FileWriter("./chicago_fairbus_k501.txt", true);
                FileWriter fileID3 = new FileWriter("./chicago_route_stops.txt", true);
                FileWriter fileID2 = new FileWriter("./chicago_busroutes_new.txt", true);
                fileID.write("geometry\\n\"{\"\"type\"\": \"\"LineString\"\", \"\"coordinates\"\": [");

                for(i=0;i<best_path.getColumnDimension();i++){
                    id = best_path.get(0, i);
                    fileID2.write((int) id);
                    fileID3.write(String.valueOf(stop_locations_data.get((int)id-1,2)) + String.valueOf(stop_locations_data.get((int)id-1, 0)));
                    fileID.write(String.valueOf(stop_locations_data.get((int)id-1,2)) + String.valueOf(stop_locations_data.get((int)id-1, 0)));
                    if(i==best_path.getColumnDimension())
                        fileID.write(",");
                }
                fileID.write("]}\"\\n");
                fileID.close();
                fileID2.close();
                fileID3.close();
                q.clear();
            }
        }
    }
}
