#include <iostream>
#include <cstdlib>
#include <cstdint>
#include <omp.h>
#include <chrono>
#include<fstream>
#include <sstream>
#include<queue>
typedef struct {
  int64_t visited;
  int64_t numneighbors;
  int64_t* neighbors;
} Node;

Node** graph;


void visit_parallel(int64_t i) {
  int64_t j, k, mark;
  #pragma omp parallel for
  for (j = 0; j < graph[i]->numneighbors; j++) {
    k = graph[i]->neighbors[j];
    #pragma omp atomic capture
    mark = graph[k]->visited++;
    if (mark == 0) {
      #pragma omp task
      visit_parallel(k);
    }
  }
}

void visit_serial(int64_t i) {
  int64_t j, k, mark;
  for (j = 0; j < graph[i]->numneighbors; j++) {
    k = graph[i]->neighbors[j];
    mark = graph[k]->visited++;
    if (mark == 0) {
      visit_serial(k);
    }
  }
}
/*
int main(int argc, char* argv[]) {
  if (argc != 3) {
    std::cerr << "Usage: ./bfs <num_nodes> <num_edges>\n";
    return 1;
  }

  int64_t num_nodes = std::atoll(argv[1]);
  int64_t num_edges = std::atoll(argv[2]);

  graph = new Node*[num_nodes];
  int64_t* degrees = new int64_t[num_nodes];
  for (int64_t i = 0; i < num_nodes; ++i) {
    degrees[i] = 0;
    graph[i] = nullptr;
  }

  int64_t* edge_from = new int64_t[num_edges];
  int64_t* edge_to = new int64_t[num_edges];

  srand(42);

  for (int64_t i = 0; i < num_edges; ++i) {
    int64_t from = rand() % num_nodes;
    int64_t to = rand() % num_nodes;
    edge_from[i] = from;
    edge_to[i] = to;
    degrees[from]++;
    degrees[to]++;
  }

  for (int64_t i = 0; i < num_nodes; ++i) {
    graph[i] = new Node;
    graph[i]->visited = 0;
    graph[i]->numneighbors = degrees[i];
    graph[i]->neighbors = new int64_t[degrees[i]];
    degrees[i] = 0;
  }

  for (int64_t i = 0; i < num_edges; ++i) {
    int64_t from = edge_from[i];
    int64_t to = edge_to[i];
    graph[from]->neighbors[degrees[from]++] = to;
    graph[to]->neighbors[degrees[to]++] = from;
  }

  delete[] degrees;
  delete[] edge_from;
  delete[] edge_to;

  auto start = std::chrono::high_resolution_clock::now();

  #pragma omp parallel
  {
    #pragma omp single
    {
      graph[0]->visited = 1;
      visit(0);
    }
  }

  #pragma omp taskwait

  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = end - start;
  std::cout << "BFS took " << elapsed.count() << " seconds.\n";

  for (int64_t i = 0; i < num_nodes; ++i) {
    delete[] graph[i]->neighbors;
    delete graph[i];
  }
  delete[] graph;

  return 0;
}*/

// Function to initialize nodes
Node* createNode() {
    Node* node = new Node;
    node->visited = 0;
    node->numneighbors = 0;
    node->neighbors = nullptr;
    return node;
}

// Function to add a neighbor to a node
void addNeighbor(Node* node, int64_t neighbor) {
    node->numneighbors++;
    int64_t* newNeighbors = new int64_t[node->numneighbors];
    for (int64_t i = 0; i < node->numneighbors - 1; i++) {
        newNeighbors[i] = node->neighbors[i];
    }
    newNeighbors[node->numneighbors - 1] = neighbor;
    delete[] node->neighbors;
    node->neighbors = newNeighbors;
}

void parallel()
{
  int numNodes=0;

  std::ifstream file("congress.edgelist");
    if (!file) {
        std::cerr << "Unable to open file" << std::endl;
        return;
    }

    // You might need to update this depending on how you determine the number of nodes.
    // For simplicity, we can initialize with a reasonable guess and resize later if needed.
    numNodes = 200; // Arbitrary number, adjust based on your graph
    graph = new Node*[numNodes];
    for (int i = 0; i < numNodes; i++) {
        graph[i] = createNode();
    }

    std::string line;
    int64_t node1, node2;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string temp;
        // Read node1, node2 and ignore the weight
        if (iss >> node1 >> node2) {
            // Resize graph if node index exceeds current size
            if (node1 >= numNodes || node2 >= numNodes) {
                int newSize = std::max(node1, node2) + 1;
                Node** newGraph = new Node*[newSize];
                for (int i = 0; i < newSize; i++) {
                    if (i < numNodes) {
                        newGraph[i] = graph[i];
                    } else {
                        newGraph[i] = createNode();
                    }
                }
                delete[] graph;
                graph = newGraph;
                numNodes = newSize;
            }
            // Add neighbors (ignore the weight part)
            addNeighbor(graph[node1], node2);
            addNeighbor(graph[node2], node1); // Assuming undirected graph
        }
    }

    file.close();

    // For debugging: print the graph
    /*for (int i = 0; i < numNodes; i++) {
        std::cout << "Node " << i << ": visited=" << graph[i]->visited
                  << ", numneighbors=" << graph[i]->numneighbors << ", neighbors={";
        for (int j = 0; j < graph[i]->numneighbors; j++) {
            std::cout << graph[i]->neighbors[j] << (j < graph[i]->numneighbors - 1 ? ", " : "");
        }
        std::cout << "}" << std::endl;
    }*/


  auto start = std::chrono::high_resolution_clock::now();
  #pragma omp parallel
  {
     #pragma omp single
     {
      graph[0]->visited = 1;
      visit_parallel(0);
    }
  }

 // #pragma omp taskwait

  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = end - start;
  std::cout << "BFS parallel took :" << elapsed.count() << " seconds.\n";
  // free the memory

  for (int64_t i = 0; i < numNodes; ++i) {
    delete[] graph[i]->neighbors;
    delete graph[i];
  }
  delete[] graph;
   // delete[] graph;

    return;
}

void serial()
{
  int numNodes=0;

  std::ifstream file("congress.edgelist");
    if (!file) {
        std::cerr << "Unable to open file" << std::endl;
        return;
    }

    // You might need to update this depending on how you determine the number of nodes.
    // For simplicity, we can initialize with a reasonable guess and resize later if needed.
    //numNodes = 200; // Arbitrary number, adjust based on your graph
    graph = new Node*[numNodes];
    for (int i = 0; i < numNodes; i++) {
        graph[i] = createNode();
    }

    std::string line;
    int64_t node1, node2;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string temp;
        // Read node1, node2 and ignore the weight
        if (iss >> node1 >> node2) {
            // Resize graph if node index exceeds current size
            if (node1 >= numNodes || node2 >= numNodes) {
                int newSize = std::max(node1, node2) + 1;
                Node** newGraph = new Node*[newSize];
                for (int i = 0; i < newSize; i++) {
                    if (i < numNodes) {
                        newGraph[i] = graph[i];
                    } else {
                        newGraph[i] = createNode();
                    }
                }
                delete[] graph;
                graph = newGraph;
                numNodes = newSize;
            }
            // Add neighbors (ignore the weight part)
            addNeighbor(graph[node1], node2);
            addNeighbor(graph[node2], node1); // Assuming undirected graph
        }
    }

    file.close();

    // For debugging: print the graph
    /*for (int i = 0; i < numNodes; i++) {
        std::cout << "Node " << i << ": visited=" << graph[i]->visited
                  << ", numneighbors=" << graph[i]->numneighbors << ", neighbors={";
        for (int j = 0; j < graph[i]->numneighbors; j++) {
            std::cout << graph[i]->neighbors[j] << (j < graph[i]->numneighbors - 1 ? ", " : "");
        }
        std::cout << "}" << std::endl;
    }*/

    // Free memory
   /* for (int i = 0; i < numNodes; i++) {
        delete[] graph[i]->neighbors;
        delete graph[i];
    }*/

  auto start = std::chrono::high_resolution_clock::now();
      graph[0]->visited = 1;
      visit_serial(0);

  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = end - start;
  std::cout << "BFS serial took   :" << elapsed.count() << " seconds.\n";
  // free the memory

  for (int64_t i = 0; i < numNodes; ++i) {
    delete[] graph[i]->neighbors;
    delete graph[i];
  }
  delete[] graph;
   // delete[] graph;

  return ;
}
int main() {
  parallel();
    serial();
    return 0;
}