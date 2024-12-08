#include<stdio.h>
#include<iostream>
#include<vector>
#include<omp.h>
#include<chrono>
#include <fstream>    // For ifstream
#include <sstream>
using namespace std;
class Graph
 {
  public:

  int V;
  vector<int> *adj;

  Graph(int ver)
  {
    V=ver;
    adj = new vector<int>[ver];
  }
  void addEdge(int u, int v)
  {
    adj[u].push_back(v);
  }
 };

void readEdgesFromFile(const string& filename, Graph &g) {
    ifstream infile(filename);
    string line;

    if (!infile) {
        cerr << "Error opening file: " << filename << endl;
        return;
    }

    while (getline(infile, line)) {
        istringstream iss(line);
        int u, v;
        string ignore;

        // Read two integers (u and v), then ignore the rest of the line
        if (iss >> u >> v) {
            g.addEdge(u, v);
        }
    }

    infile.close();
    cout << "Edges added to the graph from " << filename << endl;
}


/*

void BFS(int s, Graph &g)
{
  vector<int> edgeSize;
  vector<int> edgeOffset;
  vector<int> Edges;
  vector<int> Distance(g.V,888999888);
    printf("Size=%d\n",g.V);

  for(int i=0;i<g.V;i++)
  {
    edgeOffset.push_back(Edges.size());
    edgeSize.push_back(g.adj[i].size());
    for(int j=0;j<g.adj[i].size();j++)
    {
      Edges.push_back(g.adj[i][j]);
    }
  }


      for(int i=0; i<g.V; i++)
    {
        cout<<edgeOffset[i]<<" ";
    }
    cout<<endl;

    for(int j=0; j<Edges.size();j++)
    {
        cout<<Edges[j]<<" ";
    }
    cout<<endl;
    for(int i=0; i<g.V; i++)
    {
        cout<<edgeSize[i]<<" ";
    }
    cout<<endl;
    cout<<"**"<<endl;


    Distance[s]=0;
  int flag=1;
  int level=0;

  while(flag) {
    flag = 0;
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        if(tid < g.V && Distance[tid] == level) {
            int u = tid;
            #pragma omp parallel for
            for(int i = edgeOffset[u]; i < edgeOffset[u] + edgeSize[u]; i++) {
                int v = Edges[i];
                if(level + 1 < Distance[v]) {
                    Distance[v] = level + 1;
                    #pragma omp critical  // Ensure only one thread updates the flag at a time
                    {
                        flag = 1;
                    }
                }
            }
        }
    }
    level++;  // Increment the BFS level to continue to the next layer of nodes
}

    for(int i=0; i<g.V; i++)
    {
        cout<<"Vertex "<<i<<" : "<<Distance[i]<<endl;
    }

}

*/

void BFS_Serial(int s, Graph &g) {
    vector<int> edgeSize;
    vector<int> edgeOffset;
    vector<int> Edges;
    vector<int> Distance(g.V, 888999888);
    Distance[s] = 0;

    // Initialize edgeOffset, edgeSize, and Edges for BFS traversal
    for (int i = 0; i < g.V; i++) {
        edgeOffset.push_back(Edges.size());
        edgeSize.push_back(g.adj[i].size());
        for (int j = 0; j < g.adj[i].size(); j++) {
            Edges.push_back(g.adj[i][j]);
        }
    }

    int flag = 1;
    int level = 0;

    while (flag) {
    flag = 0;
    for (int tid = 0; tid < g.V; tid++) {
        if (Distance[tid] == level) {  // Current BFS level
            int u = tid;
            for (int i = edgeOffset[u]; i < edgeOffset[u] + edgeSize[u]; i++) {
                int v = Edges[i];
                if (level + 1 < Distance[v]) {
                    Distance[v] = level + 1;
                    flag = 1;  // Set flag to 1 if any distance is updated
                }
            }
        }
    }
    level++;  // Increment level for the next layer of BFS
}


    // Print the distances
    // for (int i = 0; i < g.V; i++) {
    //     cout << "Vertex " << i << " : " << Distance[i] << endl;
    // }
}




void BFS_Parallel(int s, Graph &g) {
    vector<int> edgeSize;
    vector<int> edgeOffset;
    vector<int> Edges;
    vector<int> Distance(g.V, 888999888);
    Distance[s] = 0;

    // Initialize edgeOffset, edgeSize, and Edges for BFS traversal
    for (int i = 0; i < g.V; i++) {
        edgeOffset.push_back(Edges.size());
        edgeSize.push_back(g.adj[i].size());
        for (int j = 0; j < g.adj[i].size(); j++) {
            Edges.push_back(g.adj[i][j]);
        }
    }

    int flag = 1;
    int level = 0;

    while (flag) {
        flag = 0;
        #pragma omp parallel for
        for (int tid = 0; tid < g.V; tid++) {
            if (Distance[tid] == level) {  // Current BFS level
                int u = tid;
                for (int i = edgeOffset[u]; i < edgeOffset[u] + edgeSize[u]; i++) {
                    int v = Edges[i];
                    if (level + 1 < Distance[v]) {
                        Distance[v] = level + 1;
                        #pragma omp atomic
                        flag |= 1;
                    }
                }
            }
        }
        level++;  // Increment level for next layer of BFS
    }

    //Print the distances
    for (int i = 0; i < g.V; i++) {
        cout << "Vertex " << i << " : " << Distance[i] << endl;
    }
}


 int main()
 {
  int numVertices = 300;  // Adjust based on the dataset
    Graph g(numVertices);
    g.addEdge(0,1);
    for(int i=0;i<150;i++)
    {
        g.addEdge(1,i);
    }
    for(int i=150;i<300;i++)
    {
        g.addEdge(2,i);
    }
    //g.addEdge(0,1);

    //readEdgesFromFile("congress.edgelist", g);

  auto Parallel_start = chrono::high_resolution_clock::now();

// Your BFS code here
BFS_Parallel(0,g);

auto Parallel_end = chrono::high_resolution_clock::now();
cout << "Parallel Time taken: " << chrono::duration_cast<chrono::nanoseconds>(Parallel_end - Parallel_start).count() << " ms" << endl;



auto start = chrono::high_resolution_clock::now();

// Your BFS code here
  BFS_Serial(1,g);

auto end = chrono::high_resolution_clock::now();
cout << "Serial Time taken: " << chrono::duration_cast<chrono::nanoseconds>(end - start).count() << " ns" << endl;

  
  
  
  cout << endl <<endl << "Serial Below " << endl;


 }