#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <stdio.h>
#include <stdlib.h>
#define netfile_path "C://Users//aabhyudai//Desktop//ibm02.net"
#define areafile_path "C://Users//aabhyudai//Desktop//ibm02.are"
#define Output_filepath "C://Users//aabhyudai//Desktop//Output02.txt"
#define bucketList_filepath "C://Users//aabhyudai//Desktop//BucketList02.txt"
#define outputNetwork_filepath "C://Users//aabhyudai//Desktop//OutputNetwork02.txt"
using namespace std;
unsigned int pins = 0, nets = 0, modules = 0, cells = 0;
unsigned int part, totalArea = 0, areaPart1 = 0, areaPart2 = 0;
unsigned int maxConnections = 0;
struct createAdjCell
{
    unsigned int cell;
    unsigned int net;
    struct createAdjCell *next;
};
struct CreateCellHeader
{
    int gain;
    bool lock_status;
    bool parti_head;
    unsigned area;
    struct createAdjCell *head; // header of adjacency list
};
struct netlist
{
    int V;
    struct CreateCellHeader *array;
};
struct bucketListCell
{
    unsigned int cell;
    struct bucketListCell *next;
};
class bucketList // For creating objects which are part of the bucketlist
{
private:
    unsigned int maxGain;
    bucketListCell *bucketArray;

public:
    bucketList(unsigned int maxConnections)
    {
        bucketArray = new bucketListCell[(2 * maxConnections) + 1];
        for (unsigned int i = (2 * maxConnections) + 1; i > 0; i--)
        {
            bucketArray[i - 1].cell = i - 1;
            bucketArray[i - 1].next = NULL;
        }
    }
    void addToList(unsigned int cell, int gain)
    {
        if (bucketArray[(maxConnections + 1 + gain)].next == NULL)
        {
            bucketListCell *newNode = new bucketListCell;
            newNode->cell = cell;
            newNode->next = NULL;
            bucketArray[(maxConnections + 1 + gain)].next = newNode;
        }
        else
        {
            struct bucketListCell *tmp;
            tmp = bucketArray[(maxConnections + 1 + gain)].next;
            while (tmp->next != NULL)
            {
                tmp = tmp->next;
            }
            bucketListCell *newNode = new bucketListCell;
            newNode->cell = cell;
            newNode->next = NULL;
            tmp->next = newNode;
        }
    }
    void outputBucketList(unsigned int maxGain)
    {
        int v;
        ofstream AdjOpFile;
        AdjOpFile.open(bucketList_filepath);
        for (v = (maxConnections + 1 - maxGain); v <= (maxConnections + 1 + maxGain); ++v)
        {
            struct bucketListCell *tmp = bucketArray[v].next;
 AdjOpFile<<"______________________________________________________________________"<<(int)(v-(
 \n";
 while(tmp!=NULL)
 {
                AdjOpFile << tmp->cell << "; ";
                tmp = tmp->next;
 }
 AdjOpFile<<"\n"<<endl;
        }
        AdjOpFile.close();
    }
 unsigned int returnNextMove(struct netlist* netlist,unsigned int maxGain,unsigned int
 {
        unsigned int leastArea, leastAreaCell;
        struct bucketListCell *tmp = bucketArray[(maxConnections + 1) + maxGain].next;
        leastArea = netlist->array[tmp->cell].area;
        leastAreaCell = tmp->cell;
        while (tmp->next != NULL)
        {
            if (leastArea > netlist->array[tmp->cell].area)
            {
                leastArea = netlist->array[tmp->cell].area;
                leastAreaCell = tmp->cell;
            }
            tmp = tmp->next;
        }
        return (leastAreaCell);
 }
 void emptyBucketList(unsigned int maxConnections)
 {
        for (unsigned int i = (2 * maxConnections) + 1; i > 0; i--)
        {
            bucketArray[i - 1].next = NULL;
        }
 }
};
struct createAdjCell *newcreateAdjCell(unsigned int cell, unsigned int net)
{
    struct createAdjCell *newNode = (struct createAdjCell *)malloc(sizeof(struct createAdjCell));
    newNode->cell = cell;
    newNode->net = net;
    newNode->next = NULL;
    return newNode;
}
struct netlist *createnetlist(unsigned int V) // Create a network based on the inputs of the ibm.net file
{
    struct netlist *netlist = (struct netlist *)malloc(sizeof(struct netlist));
    netlist->V = V;
    netlist->array = (struct CreateCellHeader *)malloc(V * sizeof(struct CreateCellHeader));
    unsigned int i;
    for (i = 0; i < V; ++i)
    {
        netlist->array[i].head = (struct createAdjCell *)newcreateAdjCell(i, 0);
        netlist->array[i].gain = 0;
        netlist->array[i].lock_status = 0;
        netlist->array[i].area = 0;
    }
    return netlist;
}
void addEdge(struct netlist *netlist, unsigned int src, unsigned int dest, unsigned int net) // joining headers
    nodes
{
    struct createAdjCell *point_srcList, *point_destList;
    point_srcList = netlist->array[src].head;
    point_destList = netlist->array[dest].head;
    struct createAdjCell *newNode = newcreateAdjCell(dest, net);
    while (point_srcList->next != NULL)
    {
        point_srcList = point_srcList->next;
    }
    point_srcList->next = newNode;
    newNode = newcreateAdjCell(src, net);
    while (point_destList->next != NULL)
    {
        point_destList = point_destList->next;
    }
    point_destList->next = newNode;
}
void writeNetlistOutput(struct netlist *netlist)
{
    unsigned int v;
    ofstream AdjOpFile;
    unsigned int counter = 0;
    AdjOpFile.open(outputNetwork_filepath);
    for (v = 0; v < netlist->V; ++v)
    {
        counter = 0;
        struct createAdjCell *tmp = netlist->array[v].head;
        AdjOpFile << "a" << v << " -> ";
        tmp = tmp->next; // Skip the first since the node is connected to itself
        while (tmp)
        {
            AdjOpFile << "a" << tmp->cell << " [" << tmp->net << "]; ";
            tmp = tmp->next;
            counter++;
        }
        AdjOpFile << "\n";
        if (maxConnections < counter)
            maxConnections = counter;
    }
    AdjOpFile.close();
}
unsigned int updateGains(struct netlist *netlist)
{
    unsigned int v, i;
    int maxGain = (-1 * maxConnections);
    bool flag = 0;
    unsigned int gain_table[maxConnections][3];
    struct createAdjCell *tmp;
    unsigned char Fs, Ts;
    for (v = 0; v <= cells; ++v) // For each cell which is not locked
    {
        Fs = 0;
        Ts = 0;
        if (netlist->array[v].lock_status == 0)
        {
            tmp = netlist->array[v].head;
            tmp = tmp->next; // Skip the first since the node is connected to itself
            unsigned char gt_row_count = 0;
            for (i = 0; i < maxConnections; i++) // Re-initializing gain table
            {
                gain_table[i][0] = 0;
                gain_table[i][1] = 0;
                gain_table[i][2] = 0;
            }
            while (tmp) // For each cell adjacent to the head cell
            {
                flag = 0;
                for (i = 0; i < gt_row_count; i++) // Scan through the gain table if the net exists
                {
                    if (gain_table[i][0] == tmp->net)
                    {
                        flag = i;
                        break;
                    }
                }
                if (flag == 0)
                {
                    gain_table[gt_row_count][0] = tmp->net;
                    // cout<<gain_table[gt_row_count][0]<<endl;
                    flag = gt_row_count;
                    gt_row_count++;
                }
                if (netlist->array[tmp->cell].parti_head == 0)
                    gain_table[flag][1]++;
                else
                    gain_table[flag][2]++;
                tmp = tmp->next;
            }
            for (i = 0; i < gt_row_count; i++)
            {
                if (netlist->array[v].parti_head == 0 && gain_table[i][1] == 0 && gain_table[i][2] > 0) // HeadCell in
                    1 and all other in partition 2 Fs++;
                else if (netlist->array[v].parti_head == 0 && gain_table[i][2] == 0 && gain_table[i][1] > 0) // HeadCell in
                    1 and all other in partition 1 Ts++;
                if (netlist->array[v].parti_head == 1 && gain_table[i][1] == 0 && gain_table[i][2] > 0) // HeadCell in
                    2 and all other in partition 2 Ts++;
                else if (netlist->array[v].parti_head == 1 && gain_table[i][2] == 0 && gain_table[i][1] > 0) // HeadCell in
                    2 and all other in partition 1 Fs++;
            }
            netlist->array[v].gain = (Fs - Ts);
            if (netlist->array[v].gain > maxGain)
                maxGain = netlist->array[v].gain;
        }
    }
    return maxGain;
}
unsigned int calculateCutSet(struct netlist *netlist)
{
    unsigned int v;
    unsigned int cutsetsize = 0;
    for (v = 0; v < netlist->V; ++v)
    {
        struct createAdjCell *tmp = netlist->array[v].head;
        tmp = tmp->next; // Skip the first since the node is connected to itself
        while (tmp)
        {
            if (netlist->array[v].parti_head != netlist->array[tmp->cell].parti_head)
                cutsetsize++;
            tmp = tmp->next;
        }
    }
    cutsetsize = cutsetsize / 2; // since each cut is counted twice (undirected list)
    return cutsetsize;
}
int main()
{
    unsigned int nRows = 0;                      // Total Number of Rows in the file
    int maxGain = 0;                             // Maximum gain in each iteration
    unsigned int rowLen = 0;                     // Number of characters in each row
    unsigned int cell_source = 0, cell_dest = 0; // to add edges
    unsigned int inNet = 0;                      // keep a count of the net through which two nodes are connected
    bool netStartsatA = 0;                       // source is a cell and not pad
    unsigned int iterations = 0;
    ofstream AdjOpFile;
    AdjOpFile.open(Output_filepath);
    std::ifstream infile(netfile_path);
    std::string line;
    //_______________________________________ .net file is read and netlist of cells is created
    while (std::getline(infile, line) && nRows <= 4)
    {
        rowLen = 0;
        unsigned char i = 0;         // counter
        while (line[rowLen] != '\0') // calculationg number of characters in the row
            rowLen++;
        char param[10];
        for (i = 0; i < rowLen; i++)
        {
            param[i] = line[i];
        }
        if (nRows == 1)
            pins = atol(param);
        else if (nRows == 2)
            nets = atol(param);
        else if (nRows == 3)
            modules = atol(param);
        else if (nRows == 4)
        {
            cells = atol(param);
            part = cells / 2;
        }
        std::istringstream iss(line);
        nRows++;
    }
    struct netlist *netlist = createnetlist(cells + 1);
    while (std::getline(infile, line))
    {
        rowLen = 0;
        unsigned char i = 0; // counter
        char nodeName[20];   // name of the pad or cell
        unsigned int cellNumber = 0;
        while (line[rowLen] != '\0') // calculationg number of characters in the row
            rowLen++;
        if (line[rowLen - 2] == '1' && line[rowLen - 3] == ' ' && line[rowLen - 4] == 's')
        {
            inNet++;
            if (line[0] == 'p')
                netStartsatA = 0;
        }
        if (line[0] == 'a')
        {
            for (i = 1; line[i] != ' '; i++)
            {
                nodeName[i - 1] = line[i];
            }
            nodeName[i - 1] = '\0';
            cellNumber = atol(nodeName);
            if (line[rowLen - 2] == '1' && line[rowLen - 3] == ' ' && line[rowLen - 4] == 's')
            {
                cell_source = cellNumber;
                netStartsatA = 1;
            }
            else
            {
                if (netStartsatA == 1)
                {
                    cell_dest = cellNumber;
                    addEdge(netlist, cell_source, cell_dest, inNet);
                }
            }
        }
        std::istringstream iss(line);
        nRows++;
    }
    //______________________________________.are file is read, areas added to corresponding
    std::ifstream infile_area(areafile_path);
    nRows = 0;
    while (std::getline(infile_area, line))
    {
        if (line[0] == 'p')
            break;
        unsigned char i = 0, j = 0; // counter
        char area[20];              // name of the pad or cell
        while (line[i] != ' ')
        {
            i++;
        }
        while (line[i] != '\0')
        {
            area[j] = line[i];
            i++;
            j++;
        }
        area[j] = '\0';
        netlist->array[nRows].area = atol(area);
        totalArea = totalArea + netlist->array[nRows].area; // total area of cells
        std::istringstream iss(line);
        nRows++;
    }
    for (unsigned int i = 0; i < cells; i++) // Initial Partitioning based on area [60-40]
    {
        areaPart1 = areaPart1 + netlist->array[i].area;
        if (areaPart1 < (0.6 * totalArea))
        {
            netlist->array[i].parti_head = 0;
            part = i;
        }
        else
        {
            netlist->array[i].parti_head = 1;
            areaPart2 = areaPart2 + netlist->array[i].area;
        }
    }
    areaPart1 = totalArea - areaPart2;
    writeNetlistOutput(netlist); // print the adjacency list representation of the above netlist
    cout << "*" << endl;
    // cout<<"Max Connections: "<<maxConnections<<endl;
    // cout<<"Number of pins: "<<pins<< endl;
    // cout<<"Number of nets: "<<nets<< endl;
    // cout<<"Number of modules: "<<modules<<endl;
    // cout<<"Number of cells: "<<cells<<endl;
    // cout<<"Total Area: "<<totalArea<<endl;
    // cout<<"Point of initial partition: between a"<<part<<" & a"<<part+1<<endl;
    // maxGain= updateGains(netlist);
    // cout<<"Initial cut-set size: "<<calculateCutSet(netlist)<<endl;
    // cout<<"Initial Partition 1 Area: "<<areaPart1<<endl;
    // cout<<"Initial Partition 2 Area: "<<areaPart2<<endl;
    // cout<<"maxGain for initial partitioning: "<<maxGain<<endl;
    AdjOpFile << AdjOpFile << AdjOpFile << "Date: 10/16/2017" << endl;
    AdjOpFile << "Project 1" << endl;
    AdjOpFile << "ESE 556, VLSI Physical and Logic Design Automation" << endl;
    AdjOpFile << "Fall 2017, Stony Brook University" << endl;
    AdjOpFile << "Author: Abhyudai Abhyudai, Chaturbhuj Rajendran, Nikolaos Maltas" << endl;
    AdjOpFile << AdjOpFile << "Number of pins: " << pins << endl;
    AdjOpFile << "Number of nets: " << nets << endl;
    AdjOpFile << "Number of modules: " << modules << endl;
    AdjOpFile << "Number of cells: " << cells << endl;
    AdjOpFile << "Total Area: " << totalArea << endl;
    AdjOpFile << "Point of initial partition: between a" << part << " & a" << part + 1 << endl;
    maxGain = updateGains(netlist);
    AdjOpFile << "Initial cut-set size: " << calculateCutSet(netlist) << endl;
    AdjOpFile << "Initial Partition 1 Area: " << areaPart1 << endl;
    AdjOpFile << "Initial Partition 2 Area: " << areaPart2 << endl;
    AdjOpFile << "maxGain for initial partitioning: " << maxGain << endl;
    AdjOpFile << bucketList bl(maxConnections);
    unsigned char flag = 0;
    // while(maxGain>0) //iterative loop for Fiduccia-Mattheyses
    while (areaPart1 > (0.58 * totalArea)) // Area Constraint
    {
        unsigned int tmp;
        flag = 0;
        iterations++;
        for (unsigned int i = 0; i < cells; i++)
        {
            if (netlist->array[i].lock_status == 0) // only if the cell is not locked, it is added to bucket list
            {
                bl.addToList(i, netlist->array[i].gain);
                flag = 1;
            }
        }
        if (flag == 0) // All nodes are locked
            break;
        // bl.outputBucketList(maxGain);
        tmp = bl.returnNextMove(netlist, maxGain, maxConnections);
        AdjOpFile << "node moved: " << tmp << endl;
        if (netlist->array[tmp].parti_head == 1)
        {
            areaPart2 = areaPart2 + netlist->array[tmp].area; // Node Moved from partition 1 to partition 2
            areaPart1 = areaPart1 - netlist->array[tmp].area;
        }
        else
        {
            areaPart1 = areaPart1 + netlist->array[tmp].area; // Node Moved from partition 2 to partition 1
            areaPart2 = areaPart2 - netlist->array[tmp].area;
        }
        netlist->array[tmp].parti_head = !(netlist->array[tmp].parti_head);
        netlist->array[tmp].lock_status = 1;
        cout << "Cut-set size: " << calculateCutSet(netlist) << endl;
        AdjOpFile << "Cut-set size: " << calculateCutSet(netlist) << endl;
        maxGain = updateGains(netlist);
        bl.emptyBucketList(maxConnections);
        cout << "maxGain after iteration " << iterations << " : " << maxGain << endl;
        AdjOpFile << "maxGain after iteration " << iterations << " : " << maxGain << endl;
        AdjOpFile << "Area of Partition 1: " << areaPart1 << endl;
        AdjOpFile << "Area of Partition 2: " << areaPart2 << endl;
        cout << AdjOpFile <<
    }
    return 0;
} // end of main