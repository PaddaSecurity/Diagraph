

#ifndef DIGRAPH_HPP
#define DIGRAPH_HPP

#include <exception>
#include <functional>
#include <list>
#include <map>
#include <utility>
#include <vector>
#include <iostream>



// DigraphExceptions are thrown from some of the member functions in the
// Digraph class template, so that exception is declared here, so it
// will be available to any code that includes this header file.

class DigraphException : public std::runtime_error
{
public:
    DigraphException(const std::string& reason);
};


inline DigraphException::DigraphException(const std::string& reason)
    : std::runtime_error{reason}
{
}



// A DigraphEdge lists a "from vertex" (the number of the vertex from which
// the edge points), a "to vertex" (the number of the vertex to which the
// edge points), and an EdgeInfo object.  Because different kinds of Digraphs
// store different kinds of edge information, DigraphEdge is a struct template.

template <typename EdgeInfo>
struct DigraphEdge
{
    int fromVertex;
    int toVertex;
    EdgeInfo einfo;
};



// A DigraphVertex includes two things: a VertexInfo object and a list of
// its outgoing edges.  Because different kinds of Digraphs store different
// kinds of vertex and edge information, DigraphVertex is a struct template.

template <typename VertexInfo, typename EdgeInfo>
struct DigraphVertex
{
    VertexInfo vinfo;
    std::list<DigraphEdge<EdgeInfo>> edges;
};



// Digraph is a class template that represents a directed graph implemented
// using adjacency lists.  It takes two type parameters:
//
// * VertexInfo, which specifies the kind of object stored for each vertex
// * EdgeInfo, which specifies the kind of object stored for each edge
//
// You'll need to implement the member functions declared here; each has a
// comment detailing how it is intended to work.
//
// Each vertex in a Digraph is identified uniquely by a "vertex number".
// Vertex numbers are not necessarily sequential and they are not necessarily
// zero- or one-based.

template <typename VertexInfo, typename EdgeInfo>
class Digraph
{
public:
    // The default constructor initializes a new, empty Digraph so that
    // contains no vertices and no edges.
    Digraph();

    // The copy constructor initializes a new Digraph to be a deep copy
    // of another one (i.e., any change to the copy will not affect the
    // original).
    Digraph(const Digraph& d);

    // The move constructor initializes a new Digraph from an expiring one.
    Digraph(Digraph&& d) noexcept;

    // The destructor deallocates any memory associated with the Digraph.
    ~Digraph() noexcept;

    // The assignment operator assigns the contents of the given Digraph
    // into "this" Digraph, with "this" Digraph becoming a separate, deep
    // copy of the contents of the given one (i.e., any change made to
    // "this" Digraph afterward will not affect the other).
    Digraph& operator=(const Digraph& d);

    // The move assignment operator assigns the contents of an expiring
    // Digraph into "this" Digraph.
    Digraph& operator=(Digraph&& d) noexcept;

    // vertices() returns a std::vector containing the vertex numbers of
    // every vertex in this Digraph.
    std::vector<int> vertices() const;

    // edges() returns a std::vector of std::pairs, in which each pair
    // contains the "from" and "to" vertex numbers of an edge in this
    // Digraph.  All edges are included in the std::vector.
    std::vector<std::pair<int, int>> edges() const;

    // This overload of edges() returns a std::vector of std::pairs, in
    // which each pair contains the "from" and "to" vertex numbers of an
    // edge in this Digraph.  Only edges outgoing from the given vertex
    // number are included in the std::vector.  If the given vertex does
    // not exist, a DigraphException is thrown instead.
    std::vector<std::pair<int, int>> edges(int vertex) const;

    // vertexInfo() returns the VertexInfo object belonging to the vertex
    // with the given vertex number.  If that vertex does not exist, a
    // DigraphException is thrown instead.
    VertexInfo vertexInfo(int vertex) const;

    // edgeInfo() returns the EdgeInfo object belonging to the edge
    // with the given "from" and "to" vertex numbers.  If either of those
    // vertices does not exist *or* if the edge does not exist, a
    // DigraphException is thrown instead.
    EdgeInfo edgeInfo(int fromVertex, int toVertex) const;

    // addVertex() adds a vertex to the Digraph with the given vertex
    // number and VertexInfo object.  If there is already a vertex in
    // the graph with the given vertex number, a DigraphException is
    // thrown instead.
    void addVertex(int vertex, const VertexInfo& vinfo);

    // addEdge() adds an edge to the Digraph pointing from the given
    // "from" vertex number to the given "to" vertex number, and
    // associates with the given EdgeInfo object with it.  If one
    // of the vertices does not exist *or* if the same edge is already
    // present in the graph, a DigraphException is thrown instead.
    void addEdge(int fromVertex, int toVertex, const EdgeInfo& einfo);

    // removeVertex() removes the vertex (and all of its incoming
    // and outgoing edges) with the given vertex number from the
    // Digraph.  If the vertex does not exist already, a DigraphException
    // is thrown instead.
    void removeVertex(int vertex);

    // removeEdge() removes the edge pointing from the given "from"
    // vertex number to the given "to" vertex number from the Digraph.
    // If either of these vertices does not exist *or* if the edge
    // is not already present in the graph, a DigraphException is
    // thrown instead.
    void removeEdge(int fromVertex, int toVertex);

    // vertexCount() returns the number of vertices in the graph.
    int vertexCount() const noexcept;

    // edgeCount() returns the total number of edges in the graph,
    // counting edges outgoing from all vertices.
    int edgeCount() const noexcept;

    // This overload of edgeCount() returns the number of edges in
    // the graph that are outgoing from the given vertex number.
    // If the given vertex does not exist, a DigraphException is
    // thrown instead.
    int edgeCount(int vertex) const;

    // isStronglyConnected() returns true if the Digraph is strongly
    // connected (i.e., every vertex is reachable from every other),
    // false otherwise.
    bool isStronglyConnected() const;

    // findShortestPaths() takes a start vertex number and a function
    // that takes an EdgeInfo object and determines an edge weight.
    // It uses Dijkstra's Shortest Path Algorithm to determine the
    // shortest paths from the start vertex to every other vertex
    // in the graph.  The result is returned as a std::map<int, int>
    // where the keys are vertex numbers and the value associated
    // with each key k is the predecessor of that vertex chosen by
    // the algorithm.  For any vertex without a predecessor (e.g.,
    // a vertex that was never reached, or the start vertex itself),
    // the value is simply a copy of the key.
    std::map<int, int> findShortestPaths(
        int startVertex,
        std::function<double(const EdgeInfo&)> edgeWeightFunc) const;

    bool checkEdge(int& fromVertex, int& toVertex) const;
    void firstThrow(int& vertex) const;
    void secondThrow(int& vertex) const;
    void thirdThrow(int& vertex) const;
    bool helper(int vertex, std::vector<int>& myVec) const;
    bool checkFind(std::vector<int>& myVec, int toCheck) const;
    bool checkSize(std::vector<int>& first, const int& second) const;

private:
    std::map<int, DigraphVertex<VertexInfo, EdgeInfo>> allThings;
};


template <typename VertexInfo, typename EdgeInfo>
Digraph<VertexInfo, EdgeInfo>::Digraph()
{
}


template <typename VertexInfo, typename EdgeInfo>
Digraph<VertexInfo, EdgeInfo>::Digraph(const Digraph& d)
{
    this->allThings = d.allThings;
}


template <typename VertexInfo, typename EdgeInfo>
Digraph<VertexInfo, EdgeInfo>::Digraph(Digraph&& d) noexcept
{
    this->allThings.swap(d.allThings);
}


template <typename VertexInfo, typename EdgeInfo>
Digraph<VertexInfo, EdgeInfo>::~Digraph() noexcept
{
    for (typename std::map<int, DigraphVertex<VertexInfo, EdgeInfo>>::iterator 
        myIter = allThings.begin();
        myIter != allThings.end(); ++myIter)
    {
        myIter->second.edges.clear();
    }
    allThings.clear();
}


template <typename VertexInfo, typename EdgeInfo>
Digraph<VertexInfo, EdgeInfo>& Digraph<VertexInfo, EdgeInfo>::operator=(const Digraph& d)
{
    this->allThings = d.allThings;
    return *this;
}


template <typename VertexInfo, typename EdgeInfo>
Digraph<VertexInfo, EdgeInfo>& Digraph<VertexInfo, EdgeInfo>::operator=(Digraph&& d) noexcept
{
    this->allThings.swap(d.allThings);
    return *this;
}


template <typename VertexInfo, typename EdgeInfo>
std::vector<int> Digraph<VertexInfo, EdgeInfo>::vertices() const
{
    std::vector<int> allVertices;
    for (typename std::map<int, DigraphVertex<VertexInfo,EdgeInfo>>::const_iterator
        myIter = allThings.begin(); myIter != allThings.end(); ++myIter)
        {
            allVertices.push_back(myIter->first);
        }
    return allVertices;
}


template <typename VertexInfo, typename EdgeInfo>
std::vector<std::pair<int, int>> Digraph<VertexInfo, EdgeInfo>::edges() const
{
    std::vector<std::pair<int, int>> allEdges;

    for (typename std::map<int, DigraphVertex<VertexInfo,EdgeInfo>>::const_iterator
        myIter = allThings.begin(); myIter != allThings.end(); ++myIter)
    {
        for (typename std::list<DigraphEdge<EdgeInfo>>::const_iterator
            iter2 = myIter->second.edges.begin(); iter2 != myIter->second.edges.end(); ++iter2)
        {
            int temp1 = iter2->fromVertex;
            int temp2 = iter2->toVertex;
            std::pair<int, int> temp3 = {temp1, temp2};
    
            allEdges.push_back(temp3);
        }
    }
    
    return allEdges;
}


template <typename VertexInfo, typename EdgeInfo>
std::vector<std::pair<int, int>> Digraph<VertexInfo, EdgeInfo>::edges(int vertex) const
{
    std::vector<std::pair<int,int>> tempEdges;

    firstThrow(vertex);

    for (typename std::list<DigraphEdge<EdgeInfo>>::const_iterator
        myIter = allThings.find(vertex)->second.edges.begin();
        myIter != allThings.find(vertex)->second.edges.end();
        ++myIter)
        {
            tempEdges.push_back(std::pair<int,int>{vertex,myIter->toVertex});
        }
    
    return tempEdges;
}


template <typename VertexInfo, typename EdgeInfo>
VertexInfo Digraph<VertexInfo, EdgeInfo>::vertexInfo(int vertex) const
{
    secondThrow(vertex);

    return allThings.at(vertex).vinfo;
}


template <typename VertexInfo, typename EdgeInfo>
EdgeInfo Digraph<VertexInfo, EdgeInfo>::edgeInfo(int fromVertex, int toVertex) const
{
    firstThrow(fromVertex);
    firstThrow(toVertex);

    if (checkEdge(fromVertex, toVertex))
    {
        std::list<DigraphEdge<EdgeInfo>> temp = allThings.at(fromVertex).edges;

        for (typename std::list<DigraphEdge<EdgeInfo>>::const_iterator myIter = 
            temp.begin(); myIter != temp.end(); ++myIter)
        {
            int temp1 = myIter->toVertex;
            if (temp1 == toVertex)
            {
                return myIter->einfo;
            }
        }
    }
    else
    {
        throw DigraphException("Edge error.");
    }
    
    return EdgeInfo{};
}


template <typename VertexInfo, typename EdgeInfo>
void Digraph<VertexInfo, EdgeInfo>::addVertex(int vertex, const VertexInfo& vinfo)
{
    thirdThrow(vertex);
    
    allThings.insert(std::pair<int, DigraphVertex<VertexInfo,EdgeInfo>>(
        vertex, DigraphVertex<VertexInfo,EdgeInfo>{vinfo}));
     
}


template <typename VertexInfo, typename EdgeInfo>
void Digraph<VertexInfo, EdgeInfo>::addEdge(int fromVertex, int toVertex, const EdgeInfo& einfo)
{
    firstThrow(fromVertex);
    firstThrow(toVertex);
   
    if (!checkEdge(fromVertex, toVertex))
    {
        std::list<DigraphEdge<EdgeInfo>>& temp1 = allThings.at(fromVertex).edges;
        temp1.push_back(DigraphEdge<EdgeInfo>{fromVertex,toVertex,einfo});
    }
    else
    {
        throw DigraphException("edge error");
    }
    
}


template <typename VertexInfo, typename EdgeInfo>
void Digraph<VertexInfo, EdgeInfo>::removeVertex(int vertex)
{
    firstThrow(vertex);
    for (typename std::map<int, DigraphVertex<VertexInfo,EdgeInfo>>::iterator
        myIter = allThings.begin(); myIter != allThings.end(); ++myIter)
    {
        auto temp2 = myIter->second.edges;
        for (typename std::list<DigraphEdge<EdgeInfo>>::iterator
            iter2 = temp2.begin(); iter2 != temp2.end(); ++iter2)
        {
            if (iter2->toVertex == vertex)
            {
                myIter->second.edges.erase(iter2);
            }
        }
    }
    allThings.erase(vertex);
}


template <typename VertexInfo, typename EdgeInfo>
void Digraph<VertexInfo, EdgeInfo>::removeEdge(int fromVertex, int toVertex)
{
    firstThrow(fromVertex);
    firstThrow(toVertex);

    if (checkEdge(fromVertex, toVertex))
    {

        for (typename std::list<DigraphEdge<EdgeInfo>>::iterator myIter = allThings.at(fromVertex).edges.begin(); 
            myIter != allThings.at(fromVertex).edges.end(); ++myIter)
        {
            if (myIter->toVertex == toVertex)
            {
                // std::cout << "(" << myIter->fromVertex << ", " << myIter->toVertex << ")" << std::endl;
                allThings.at(fromVertex).edges.erase(myIter);
                break;
            }
        }       
    }
    else
    {
        throw DigraphException("Edge does not exist");
    }
}


template <typename VertexInfo, typename EdgeInfo>
int Digraph<VertexInfo, EdgeInfo>::vertexCount() const noexcept
{
    return allThings.size();
}


template <typename VertexInfo, typename EdgeInfo>
int Digraph<VertexInfo, EdgeInfo>::edgeCount() const noexcept
{
    if (allThings.size() == 0)
    {
        return 0;
    }
    unsigned int totalEdges = 0;

    for (typename std::map<int, DigraphVertex<VertexInfo, EdgeInfo>>::const_iterator
        myIter = allThings.begin(); myIter != allThings.end(); ++myIter)
    {
        auto temp1 = myIter->second.edges;
        totalEdges += temp1.size();
    }

    return totalEdges;
}


template <typename VertexInfo, typename EdgeInfo>
int Digraph<VertexInfo, EdgeInfo>::edgeCount(int vertex) const
{
    int returnVal = 0;
    if (allThings.count(vertex) != 1)
    {
        throw DigraphException("Vertex not found in map");
    }
    returnVal = allThings.at(vertex).edges.size();
    return returnVal;
}


template <typename VertexInfo, typename EdgeInfo>
bool Digraph<VertexInfo, EdgeInfo>::isStronglyConnected() const
{       
    for (typename std::map<int, DigraphVertex<VertexInfo,EdgeInfo>>::const_iterator
        myIter = allThings.begin(); myIter != allThings.end(); ++myIter)
    {
        std::vector<int> temp;
        helper(myIter->first, temp);
        
        for (unsigned i = 0; i < temp.size(); ++i)
        {
            std::cout << temp[i] << " ";
        }
        std::cout << std::endl;

        if (temp.size() != allThings.size())
        {
            return false;
        }
    }

    return true; 
}

//checkCon helpers
template<typename VertexInfo, typename EdgeInfo>
bool Digraph<VertexInfo,EdgeInfo>::helper(int vertex, std::vector<int>& allLeft) const
{
    for (typename std::list<DigraphEdge<EdgeInfo>>::const_iterator
        myIter = allThings.at(vertex).edges.begin(); 
        myIter != allThings.at(vertex).edges.end();
        ++myIter)
    {
        if (std::find(allLeft.begin(), allLeft.end(), myIter->toVertex) == allLeft.end())
        {
            allLeft.push_back(myIter->toVertex);
        }
        else if (allLeft.size() < allThings.size() - 1)
        {
            helper(myIter->toVertex, allLeft);
        }
        else
        {
            return true;
        }
    }

    return false;
}
template<typename VertexInfo, typename EdgeInfo>
bool Digraph<VertexInfo,EdgeInfo>::checkFind(std::vector<int>& myVec, int toCheck) const
{
    if (std::find(myVec.begin(), myVec.end(), toCheck) != myVec.end())
    {
        return true;
    }
    return false;
}
template<typename VertexInfo, typename EdgeInfo>
bool Digraph<VertexInfo,EdgeInfo>::checkSize(std::vector<int>& first, const int& second) const
{
    return first.size() < second - 1 ? true : false;
}

template <typename VertexInfo, typename EdgeInfo>
std::map<int, int> Digraph<VertexInfo, EdgeInfo>::findShortestPaths(
    int startVertex,
    std::function<double(const EdgeInfo&)> edgeWeightFunc) const
{

    std::map<int, int> xxx;
    std::map<int, double> shortestDistance;

    for (typename std::map<int, DigraphVertex<VertexInfo,EdgeInfo>>::const_iterator
        myIter = allThings.begin(); myIter != allThings.end(); ++myIter)
    {
        int temp = myIter->first;
        if (temp == startVertex)
        {
            shortestDistance[startVertex] = 0;
        }
        else
        {
            shortestDistance[temp] = INT_MAX;
        }
    }  

   
    xxx[startVertex] = startVertex;
    std::vector<int> allUnvisited = vertices();

    while (allUnvisited.size() != 0)
    {

        int temp2 = 0;
        double m = INT_MAX;
        int mm = 0;
    
        for(unsigned i = 0; i < allUnvisited.size(); i++)
        {
            int bob2 = allUnvisited[i];
            double bob = shortestDistance[bob2];
            
            if (bob < m)
            {
                m = bob;
                mm = bob2;
            }  
        }

        temp2 = mm;
        try
        {
            std::vector<int>::iterator it = std::find(allUnvisited.begin(), allUnvisited.end(), temp2);
            allUnvisited.erase(it);
        }
        catch(...)
        {
            std::cout << "oops" << std::endl;
        }
        
        double distance = 0.0;
        double temp3 = 0.0;

        for(typename std::list<DigraphEdge<EdgeInfo>>::const_iterator 
            myIter = allThings.at(temp2).edges.begin();
            myIter != allThings.at(temp2).edges.end();
            ++myIter)
        {

            const int& temp4 = myIter->toVertex;
            distance = shortestDistance[temp2];
            try
            {
                temp3 = edgeWeightFunc(myIter->einfo);
            }
            catch(...)
            {
                std::cout << "function throwing error" << std::endl;
            }
            
            distance += temp3;

            double& finalDist = shortestDistance[temp4];

            if(distance < finalDist)
            {
                finalDist = distance;
                xxx[temp4] = temp2;
            }
        }
    }

    return xxx;   
}

template <typename VertexInfo, typename EdgeInfo>
bool Digraph<VertexInfo,EdgeInfo>::checkEdge(int& fromVertex, int& toVertex) const
{
    std::vector<std::pair<int,int>> tempVec = edges(fromVertex);
    std::pair<int,int> myVal = std::pair<int,int>{fromVertex,toVertex};

    if (std::find(tempVec.begin(), tempVec.end(), myVal) != tempVec.end())
    {
        // std::cout << "***" << std::endl;
        // std::cout << "(" << myVal.first << ", " << myVal.second << ")" << std::endl;
        // std::cout << "***" << std::endl;

        return true;
    }
    return false;
}

template<typename VertexInfo, typename EdgeInfo>
void Digraph<VertexInfo,EdgeInfo>::firstThrow(int& vertex) const
{
    if (allThings.count(vertex) != 1)
    {
        throw DigraphException("vertex error");
    }
}
template<typename VertexInfo, typename EdgeInfo>
void Digraph<VertexInfo,EdgeInfo>::secondThrow(int& vertex) const
{
    if (allThings.count(vertex) == 0)
    {
        throw DigraphException("vertex error");
    }
}
template<typename VertexInfo, typename EdgeInfo>
void Digraph<VertexInfo,EdgeInfo>::thirdThrow(int& vertex) const
{
    if (allThings.count(vertex) != 0)
    {
        throw DigraphException("vertex error");
    }
}



#endif

