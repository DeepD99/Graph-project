/*main.cpp*/

//
// Deep Diganvker
// U. of Illinois, Chicago
// CS 251: Spring 2020
// Project #07: open street maps, graphs, and Dijkstra's alg
// 
// If you are reading this, please write me a cs joke :) 
// 
// References:
// TinyXML: https://github.com/leethomason/tinyxml2
// OpenStreetMap: https://www.openstreetmap.org
// OpenStreetMap docs:  
//   https://wiki.openstreetmap.org/wiki/Main_Page
//   https://wiki.openstreetmap.org/wiki/Map_Features
//   https://wiki.openstreetmap.org/wiki/Node
//   https://wiki.openstreetmap.org/wiki/Way
//   https://wiki.openstreetmap.org/wiki/Relation
//

#include <iostream>
#include <iomanip>  /*setprecision*/
#include <string>
#include <vector>
#include <map>
#include <set>
#include <queue>
#include <utility>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <limits>
#include <stack>

#include "graph.h"
#include "dist.h"
#include "osm.h"

using namespace std;
using namespace tinyxml2;

// Declare INF
double INF = numeric_limits<double>::max();

// 
// fileStats
// 
// Helper function that is used to print the data of the map
// Nodes, Footways, and Buildings
// 
void fileStats(graph<long long, double>& G, map<long long, Coordinates>& Nodes, 
               vector<FootwayInfo>& Footways, vector<BuildingInfo>& Buildings, string def_filename, string filename)
{    
  // Beginning of the program
  cout << "** Navigating UIC open street map **" << endl;
  cout << endl;
  cout << std::setprecision(8);

  // Declares map.osm
  def_filename = "map.osm";

  // Prompts the user to open map.osm
  cout << "Enter map filename> ";
  getline(cin, filename);

  // If the name of file matches, then proceed to the program
  if (filename == "")
    filename = def_filename;
    
  // Declare an instance of class XMLDocument:
  XMLDocument xmldoc; 
    
  // Load XML-based map file
  if (!LoadOpenStreetMap(filename, xmldoc))
  {
    cout << "**Error: unable to load open street map." << endl;
    cout << endl;
    return;
  }  
  
  // Read the nodes, which are the various known positions on the map:
  unsigned nodeCount = ReadMapNodes(xmldoc, Nodes);

  // Read the footways, which are the walking paths:
  unsigned footwayCount = ReadFootways(xmldoc, Footways);
    
  // Read the university buildings:
  unsigned buildingCount = ReadUniversityBuildings(xmldoc, Nodes, Buildings);
  
  // Test:
  assert(nodeCount == Nodes.size());
  assert(footwayCount == Footways.size());
  assert(buildingCount == Buildings.size());  
    
  // Prints the data of the map:
  cout << endl;
  cout << "# of nodes: " << Nodes.size() << endl;
  cout << "# of footways: " << Footways.size() << endl;
  cout << "# of buildings: " << Buildings.size() << endl;

  // Iterate through the Nodes map in order to add the vertices to the graph
  for(auto addV : Nodes)
      G.addVertex(addV.first);
    
  // Iterate through the footways vector:
  for(auto footWay : Footways)
  {
      for(size_t footWayNode = 0; footWayNode < footWay.Nodes.size()-1; footWayNode++)
      {
          // Call the distBetween2Points function that finds the distance between the departure and the arrival
          double travelDistance = distBetween2Points(Nodes[footWay.Nodes[footWayNode]].Lat, 
                                                     Nodes[footWay.Nodes[footWayNode]].Lon,
                                                     Nodes[footWay.Nodes[footWayNode+1]].Lat, 
                                                     Nodes[footWay.Nodes[footWayNode+1]].Lon);
          // Add the edges appropriately in the graph to compute Dijkstra's algorithm
          G.addEdge(footWay.Nodes[footWayNode], footWay.Nodes[footWayNode+1], travelDistance);
          G.addEdge(footWay.Nodes[footWayNode+1], footWay.Nodes[footWayNode], travelDistance);  
      }
  }
  
  // Output the vertices and edges of the map:
  cout << "# of vertices: " << G.NumVertices() << endl;
  cout << "# of edges: " << G.NumEdges() << endl;
  cout << endl;
}

//
// startPoint
// 
// Check for start building in main when user inputs an invalid input
// Gives the correct start building
// Gives the correct coordinates for the start building
// 
bool startPoint(double iterateStart, double& startLat, double& startLon, vector<BuildingInfo> Buildings, string startBuilding)
{
    // Go through all the buildings for start:
    for(iterateStart = 0; iterateStart < Buildings.size(); iterateStart++)
    {
       // Obtain the first string of the userInput
       auto startFullName      = Buildings.at(iterateStart).Fullname;
       // First string of the input
       auto storeStartFullName = startFullName.substr(0, Buildings.at(iterateStart).Fullname.find(' '));
       // Obtain the abbrev from the input, if given
       auto startAbbrev        = startBuilding.substr(0, startBuilding.find(' '));
        
       // Compare the fullname with startAbbrev
       // Check if the startBuilding is equal to the abbrev --> Getting the correct initial position
       if(storeStartFullName.compare(startAbbrev) == 0 || startBuilding == Buildings.at(iterateStart).Abbrev)
       {
           // Prompt
           cout << "Starting point: " << endl;
           
           // Store the latitudes and longitudes of the start distance:
           startLat = Buildings.at(iterateStart).Coords.Lat;
           startLon = Buildings.at(iterateStart).Coords.Lon;
           
           // Prints the coordinates:
           cout << " " << Buildings.at(iterateStart).Fullname << endl;
           cout << " " << "(" << startLat << ", " << startLon << ")" << endl;
           return true;
       }     
    }
            
    // If not starting building found, then print statement:
    cout << "Start building not found" << endl;
    return false;
}

//
// destPoint
// 
// Gives the correct destination building
// Gives the correct coordinates
// 
bool destPoint(double iterateDest, double& endLat, double& endLon, vector<BuildingInfo> Buildings, bool& checkStart, string destBuilding)
{
    // Go through all the buildings for end:
    for(iterateDest = 0; iterateDest < Buildings.size(); iterateDest++)
    {
       // Obtain the first string of the userInput
       auto destFullName      = Buildings.at(iterateDest).Fullname;
       // First string of the input
       auto storeDestFullName = destFullName.substr(0, Buildings.at(iterateDest).Fullname.find(' '));
       // Obtain the abbrev from the input, if given
       auto destAbbrev        = destBuilding.substr(0, destBuilding.find(' '));
        
       // Compare the fullname with startAbbrev
       // Check if the destBuilding is equal to the abbrev --> Getting the correct final position
       if(storeDestFullName.compare(destAbbrev) == 0 || destBuilding == Buildings.at(iterateDest).Abbrev)
       {
           // Prompt
           cout << "Destination point: " << endl;
           
           // Store the latitudes and longitudes of the end distance:
           endLat = Buildings.at(iterateDest).Coords.Lat;
           endLon = Buildings.at(iterateDest).Coords.Lon;
           
           // Prints the coordinates:
           cout << " " << Buildings.at(iterateDest).Fullname << endl;
           cout << " " << "(" << endLat << ", " << endLon << ")" << endl;
           // Spacing:
           if(checkStart)
               cout << endl;
           // Return true, since we have a valid, destination building:
           return true;
       }
    }   
	cout << "Destination building not found" << endl;
    return false;
}


//
// Ordering
// 
// Used for the bool function that handles the sorting algorithm
// for the priority queue
//
class ordering
{
   public:
      bool operator() (const pair<long long, double>& from, const pair<long long, double>& to)
      {
         // If you don't want to go all the way to the end of the queue, then return false
         if(from.second < to.second)
            return false;
         // If you want to go all the way to the end of the queue, then return true
         else if(from.second > to.second)
            return true;
         // Else, then return the largest, front-most key, in the queue
         else
            return from.first > to.first;
      }
};

//
// helperForCheckVisited
// 
// Returns true if there are no infinities in the sort
// Returns false if there are infinities in the sort
// 
bool helperForCheckVisited (long long nameTheVertex, vector<long long> storedVertices)
{
   // Iterate through the vector
   for(size_t iterateVector = 0; iterateVector < storedVertices.size(); iterateVector++)
   {
      // If it matches the name, then return true:
      if(storedVertices[iterateVector] == nameTheVertex)
         return true;
   }
   
   // Else, return false:
   return false;
}

//
// Dijkstra
//
// Performs Dijkstra's shortest weighted path algorithm from
// the given start vertex. 
// Returns a map of predecessors from the edge
// Prints the shortest path between two buildings
//
map<long long, long long> Dijkstra(graph<long long, double>& G, long long& startV, long long& endV,
                                   map<long long, double>& distances, 
                                   map<long long, long long> tracePathForPredecessors,
                                   stack<long long>& printShortestPath, double& saveDistance)
{
  // Every visited vertex is stored here
  vector<long long> visited;
  // What need to be visited
  priority_queue <pair<long long, double>, 
                 vector<pair<long long, double>>,
                 ordering> unvisitedQueue;
  // Calls the vertices for further operations:
  vector<long long> vertices = G.getVertices();
  // Set that makes sure that all vertices are visited
  set<long long>    checkVertices;
  // Set that get adjacent vertices
  set<long long>    adjacentVertices;
      
  // Loop through all the vertices  
  for(long long curV : vertices)
  {
     // Set all the values of the vertices to be INF
     distances[curV] = INF;
     // Initialize the predecessor to 0
     tracePathForPredecessors[curV] = 0;  
     // Push it into the unvisited queue in order to visit them step by step in the trace
     unvisitedQueue.push(make_pair(curV, INF));
  }
   
  // startV has a distance of 0 from itself and for predecessor:
  distances[startV] = 0;  
  unvisitedQueue.push(make_pair(startV, 0));
  
  // Iterate until the Priority Queue is empty:
  while (!unvisitedQueue.empty())
  {
     // Visit vertex with minimum distance from startV
     pair<long long, double> currentVertex = unvisitedQueue.top();
     unvisitedQueue.pop();
     
     // Calls the function that checks if a vertex is visited
     bool helperVisited = helperForCheckVisited(currentVertex.first, visited);
     
     // Check if there are any infinites in the sort
     if((currentVertex.second) == INF)
        break;
    
     // Check if there are any duplicates in the sort
     else if(helperVisited)
        continue;
      
     // Else, then insert the vertex to the set
     else
        checkVertices.insert(currentVertex.first);
      
     // Call the neighbors function to iterate through the set
     adjacentVertices = G.neighbors(currentVertex.first);
           
     // Iterate until you get all the neighbors
     for(long long adjacentVertex : adjacentVertices)
     {
        // Declare edgeWeight:
        double edgeWeight = 0.0;
        // Get the weight of the edge
        G.getWeight(currentVertex.first, adjacentVertex, edgeWeight);
        // Add it using the map and the weight of the edge to compute the alternativePathDistance:
        double alternativePathDistance = distances[currentVertex.first] + edgeWeight;
        
        // If shorter path from startV to adjV is found
        // Update adjV's distance and predecessor
        // Save the distance between two buildings
        // Store predecessor
        if(alternativePathDistance < distances[adjacentVertex])
        {
           distances[adjacentVertex] = alternativePathDistance;
           tracePathForPredecessors[adjacentVertex] = currentVertex.first;
           saveDistance = (distances[endV] - distances[startV]); 
           unvisitedQueue.push(make_pair(adjacentVertex, alternativePathDistance));
        }
     }
  }   
  
  // Push into the stack:
  printShortestPath.push(endV);  
    
  // Pushing the values to the tracePathForPredecessors map:
  // Sent to nearestStartAndDest in order to print the shortest path using
  // Dijkstra's algorithm
  while(tracePathForPredecessors[endV])
  {
     // Set the endV to the specific index of tracePathForPredecessors:
     endV = tracePathForPredecessors[endV];
     // Push into the stack:
     printShortestPath.push(endV);
  }

  // Return the map of the predecessor:
  return tracePathForPredecessors;
}

//
// nearestStartAndDest
// 
// Gets the nearest starting point and ending point of a building
// Gives the ID from the Coordinates struct
// Gives the coordinates for the Coordinates struct
// 
void nearestStartAndDest(vector<FootwayInfo> Footways, map<long long, Coordinates> Nodes,
                         double endOfNodes, double getDistanceForStart, double getDistanceForEnd,
                         double& startLat, double& startLon, double& endLat, double& endLon, 
                         double& initialLat, double& initialLon, double& finalLat, double& finalLon,
                         double endOfAllNodes, double getAllDistanceForStart, double getAllDistanceForEnd,
                         long long& initialID, long long& finalID, graph<long long, double>& G, 
                         stack<long long>& printShortestPath, double& saveDistance)
{        
    // Iterate through Footways:
    // For each of those, check the nodes vector in the struct FootwayInfo:
    // Gets the first index of the distance function:
    for(auto iterateFootForBuildings : Footways)
    {
        // Iterate until the end of the Nodes vector:
        for(endOfNodes = 0; endOfNodes < iterateFootForBuildings.Nodes.size(); endOfNodes++)
        {
            // Get the starting distance using the distBetween2Points function:
            getDistanceForStart = distBetween2Points(startLat, 
                                                     startLon, 
                                                     Nodes.at(iterateFootForBuildings.Nodes.at(0)).Lat,
                                                     Nodes.at(iterateFootForBuildings.Nodes.at(0)).Lon);
            
            // Get the ending distance using the distBetween2Points function:
            getDistanceForEnd   = distBetween2Points(endLat,
                                                     endLon,
                                                     Nodes.at(iterateFootForBuildings.Nodes.at(0)).Lat,
                                                     Nodes.at(iterateFootForBuildings.Nodes.at(0)).Lon);  
        }   
    }
     
    // Iterate through Footways:
    // Gets the first index from the previous for loop in order to compare for the previous location
    // Gets the shortest/nearest starting and ending node:
    for(auto iterateAllFootForBuildings : Footways)
    {
        // Iterate until the end of the nodes vector:
        for(endOfAllNodes = 0; endOfAllNodes < iterateAllFootForBuildings.Nodes.size(); endOfAllNodes++)
        {
            // Get the starting distance using the distBetween2Points function:
            getAllDistanceForStart = distBetween2Points(startLat, 
                                                        startLon, 
                                                        Nodes.at(iterateAllFootForBuildings.Nodes.at(endOfAllNodes)).Lat,
                                                        Nodes.at(iterateAllFootForBuildings.Nodes.at(endOfAllNodes)).Lon);
           
            // Get the ending distance using the distBetween2Points function:
            getAllDistanceForEnd   = distBetween2Points(endLat,
                                                        endLon,
                                                        Nodes.at(iterateAllFootForBuildings.Nodes.at(endOfAllNodes)).Lat,
                                                        Nodes.at(iterateAllFootForBuildings.Nodes.at(endOfAllNodes)).Lon);   
           
            // If the starting distance from the first index is less than the maximum distance, then update the node:
            if(getDistanceForStart > getAllDistanceForStart)
            {
                // Update the node
                getDistanceForStart = getAllDistanceForStart;
               
                // Update the data values: ID, Lat, and Lon:
                initialID  = Nodes.at(iterateAllFootForBuildings.Nodes.at(endOfAllNodes)).ID;
                initialLat = Nodes.at(iterateAllFootForBuildings.Nodes.at(endOfAllNodes)).Lat;
                initialLon = Nodes.at(iterateAllFootForBuildings.Nodes.at(endOfAllNodes)).Lon;
            }
                      
            // If the ending distance from the first index is less than the maximum distance, then update the node:
            if(getDistanceForEnd > getAllDistanceForEnd)
            {
                // Update the node
                getDistanceForEnd = getAllDistanceForEnd;
               
                // Update the data values: ID, Lat, and Lon:
                finalID  = Nodes.at(iterateAllFootForBuildings.Nodes.at(endOfAllNodes)).ID;
                finalLat = Nodes.at(iterateAllFootForBuildings.Nodes.at(endOfAllNodes)).Lat;
                finalLon = Nodes.at(iterateAllFootForBuildings.Nodes.at(endOfAllNodes)).Lon;
            }
        }
    }
   
    // Printing the nearest starting node accordingly:
    cout << "Nearest start node: " << endl;
    cout << " " << initialID << endl;
    cout << " " << "(" << initialLat << ", " << initialLon << ")" << endl;
      
    // Printing the nearest destination node accordingly:
    cout << "Nearest destination node: " << endl;
    cout << " " << finalID << endl;
    cout << " " << "(" << finalLat   << ", " << finalLon   << ")" << endl;
    cout << endl;
      
    // Starting the Dijkstra's algorithm
    cout << "Navigating with Dijkstra..." << endl;

    // Declaring the maps to call the dijkstra function:
    map<long long, long long> storedPredecessors;
    map<long long,    double> receiveDistance;  

    // Check if both buildings are valid:
    if(initialID != finalID)
    {
        // Call the Dijkstra function
        Dijkstra(G, initialID, finalID, receiveDistance, storedPredecessors, printShortestPath, saveDistance);

        // Until the map of the recieve distance is not infinity:
        if(receiveDistance[finalID] != INF)
        {  
           // Distance between two buildings:
           cout << "Distance to dest: " << saveDistance << " miles" << endl;

           // Print the path:
           cout << "Path: ";   
           while(!printShortestPath.empty())
           {
               // Iterating until we get the full path from the predecessor map:
               // Find the top element of the stack:
               long long printNode = printShortestPath.top();     
               // Print the top element:
               cout << printNode;
               // Pop the element:
               printShortestPath.pop();
               // Until the path is empty, then keep having arrows point to another node
               if(printShortestPath.empty() == storedPredecessors[finalID])
                   cout << "->";
           }       
           // Spacing:
           cout << endl;
        }
        // If an invalid building is inputted, then return an error message:
        // Sorry, destination unreachable
        else
           cout << "Sorry, destination unreachable" << endl;
    }
    // If the buildings are the same:
    else
    {
       // Print the destination distance, which will be zero
       cout << "Distance to dest: " << saveDistance << " miles" << endl;
       // Print the path:
       cout << "Path: ";   
       cout << initialID << endl;
    }   
}

//////////////////////////////////////////////////////////////////
//
// main
//
int main()
{
  map<long long, Coordinates>  Nodes;                // maps a Node ID to it's coordinates (lat, lon)
  vector<FootwayInfo>          Footways;             // info about each footway, in no particular order
  vector<BuildingInfo>         Buildings;            // info about each building, in no particular order
  graph<long long, double>     G;                    // graph that is used to add the edge in a map
  string                       def_filename;         // used to get map.osm
  string                       filename;             // asked from the user for a map file
  stack<long long>             printShortestPath;    // store the path of a trip from one building to another
	
    
  // Call the stats function that tells the data of the map
  // Foundation of the program starts here
  fileStats(G, Nodes, Footways, Buildings, def_filename, filename);
    
  // Navigation from building to building
  string startBuilding, destBuilding;

  // Prompting the user for a start building:
  cout << "Enter start (partial name or abbreviation), or #> ";
  getline(cin, startBuilding);

  // While the start building is not equal to '#':
  while (startBuilding != "#")
  {
    // Prompting the user for their destination building:
    cout << "Enter destination (partial name or abbreviation)> ";
    getline(cin, destBuilding);    
      
    // Create four integers to be used by the for-loops
    size_t iterateStart, iterateDest, endOfNodes, endOfAllNodes;
    iterateStart = iterateDest = endOfNodes = endOfAllNodes = 0;
    
    // Create four double variables to store the latitude and longitude of 
    // the start and end distance of the buildings:
    double startLat, startLon, endLat, endLon;
    startLat = startLon = endLat = endLon = 0.0;
      
    // Create two double variables for calling the distance function 
    // - Iterating for one node in Footways
    double getDistanceForStart, getDistanceForEnd;
    getDistanceForStart = getDistanceForEnd = 0.0;
      
    // Create two double variables for calling the distance function 
    // - Iterating for all nodes in Footways
    double getAllDistanceForStart, getAllDistanceForEnd;
    getAllDistanceForStart = getAllDistanceForEnd = 0.0;
      
    // Create two variables that stores the ID of the start and end 
    // once it is updated by the nearest distance:
    long long initialID, finalID;
    initialID = finalID = 0;
    
    // Create four variables that represents the LAT and LON for start and end 
    // after it is being updated by the nearest distance:
    double initialLat, initialLon, finalLat, finalLon;
    initialLat = initialLon = finalLat = finalLon = 0.0;
    
    // Store the distance between two buildings:
    double saveDistance = 0.0;
      
    // Call the helper function that prints out the starting and ending points with coordinates
    // Used to check if there is a start building:
    bool checkStart = startPoint(iterateStart, startLat, startLon, Buildings, startBuilding);
    
    // Used to check if there is a destination building:
    bool checkDest  = destPoint(iterateDest, endLat, endLon, Buildings, checkStart, destBuilding); 
              
      
    if(checkStart)
    {
        if(checkDest)
        {
            // Call the helper function named nearestStartAndDest that:
            // prints the nearest starting and destination node:
            // prints the path using the dijkstra algorithm
            // 
            nearestStartAndDest(Footways, Nodes, 
                                endOfNodes, getDistanceForStart, getDistanceForEnd,
                                startLat, startLon, endLat, endLon, 
                                initialLat, initialLon, finalLat, finalLon,
                                endOfAllNodes, getAllDistanceForStart, getAllDistanceForEnd,
                                initialID, finalID, G, printShortestPath, saveDistance);
        }
    }
      
    // another navigation?
    cout << endl;
    cout << "Enter start (partial name or abbreviation), or #> ";
    getline(cin, startBuilding);  
  }

  cout << "** Done **" << endl;

  return 0;
}
