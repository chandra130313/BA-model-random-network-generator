#include<iostream>
#include<string>
#include<map>
#include<list>
#include<stack>
#include<vector>
#include <sys/time.h>
#include <cstdlib>
#include <fstream>
#include<algorithm>
#include<limits.h>
#define RANDOM_FAILURE 1
#define DET_ATTACK 2
using namespace std;
class BasicFunction
{
    public:
    void InitRand();
    int RandNum();
};
void BasicFunction::InitRand()
{
    struct timeval tm;
    gettimeofday(&tm , NULL);
    srand(tm.tv_sec * 1000000 + tm.tv_usec);
}

int BasicFunction::RandNum()
{
    int temp;
    struct timeval tm;

    gettimeofday(&tm, NULL);
    temp = 1000000 * (rand()/(RAND_MAX + 1.0));

    return temp;
}
class Attack
{
    public:
    int apply_attack;
};
class BANetwork
{
    private:
    int m_nEdges;
    int m_nMaxDeg;
    float m_fEpsilon;

    float **m_pCorrMatrix;
    float m_fAssortCoeff;
    BasicFunction m_basfunc;
    vector<int> m_vNodeDegree;
    vector<float> m_vDegDist;

    void AddNewNode(int nNodeID, int nDegreeArrivingNodes);
    bool CanEdgeExist(int nNode1, int nNode2);

    public:
    BANetwork(int nNodes, int nInitNetworkSize, int nDegreeArrivingNodes, float fEpsilon);
    ~BANetwork();
    multimap<int ,int> m_mmEdgeList;
    int m_nNodes;
    float GetAssortCoeff();
    vector<float>& GetDegDist();
    void degree_distribution(char *fname);
    void clustering_coefficient();
    int   FindMaxDeg();
    int combination(int m,int n);
    void diameter();
    int factorial(int x);
    void   GenerateCorrelationMatrix();
    float  AssortativityCoeff();

    void   GenerateDegDist();

    void    ApplyAttack(int choice, float del_fr);
    void ApplyAttack_deterministically(int choice,float del_fr);

    void   WriteEdgeListtoFile(char *fname);
    void   WriteNodeDegreetoFile(char *fname);
    void   WriteDegDisttoFile(char *fname);
};

class Graph
{
    int V;
    list<int> adj[1000];

    void fillOrder(int v, bool visited[], stack<int> &Stack);

     void DFSUtil(int v, bool visited[]);
public:
    Graph();
    vector<int> output;
    int count;
    int maximum;
    void addEdge(int v, int w);
    void printSCCs(BANetwork &BANet,float fraction);

    Graph getTranspose(BANetwork &BANet);
};
Graph::Graph()
{
   count=0;
   maximum=INT_MIN;
}
void Graph::DFSUtil(int v, bool visited[])
{

    visited[v] = true;
    count++;
    list<int>::iterator i;
    for (i = adj[v].begin(); i != adj[v].end(); ++i)
        if (!visited[*i])
          DFSUtil(*i, visited);

}

Graph Graph::getTranspose(BANetwork &BANet)
{
    Graph g;
    multimap<int,int>::iterator v;
    for (v=BANet.m_mmEdgeList.begin();v!=BANet.m_mmEdgeList.end();v++)
    {

        list<int>::iterator i;
        for(i = adj[(*v).first].begin(); i != adj[(*v).first].end(); ++i)
        {
            g.adj[*i].push_back((*v).first);
        }
    }
    return g;
}

void Graph::addEdge(int v, int w)
{
    adj[v].push_back(w);
}

void Graph::fillOrder(int v, bool visited[], stack<int> &Stack)
{

    visited[v] = true;

    list<int>::iterator i;
    for(i = adj[v].begin(); i != adj[v].end(); ++i)
        if(!visited[*i])
            fillOrder(*i, visited, Stack);

        Stack.push(v);
}


 void Graph::printSCCs(BANetwork &BANet,float fraction)
{
    stack<int> Stack;

    bool *visited = new bool[V];


    for (multimap<int,int>::iterator v=BANet.m_mmEdgeList.begin();v!=BANet.m_mmEdgeList.end();v++)
        visited[(*v).first] = false;

        for (multimap<int,int>::iterator v=BANet.m_mmEdgeList.begin();v!=BANet.m_mmEdgeList.end();v++)
        if(visited[(*v).first] == false)
            fillOrder((*v).first, visited, Stack);

        Graph gr = getTranspose(BANet);


    for (multimap<int,int>::iterator v=BANet.m_mmEdgeList.begin();v!=BANet.m_mmEdgeList.end();v++)
        visited[(*v).first] = false;


    while (Stack.empty() == false)
    {

        int v = Stack.top();
        Stack.pop();

        vector<int> output;
        if (visited[v] == false)
        {
           gr.DFSUtil(v, visited);
           if(gr.maximum<gr.count)
           {
               gr.maximum=gr.count;
               gr.count=0;
           }

        }
    }
  fstream fwr;
  fwr.open("deterministic_connect.txt",ios::out|ios::app);
  cout<<gr.maximum;
  fwr<<(float)gr.maximum/(float)BANet.m_nNodes<<" "<<fraction<<endl;

}

BANetwork::BANetwork(int nNodes, int nInitNetworkSize, int nDegreeArrivingNodes, float fEpsilon=0)
{
    m_basfunc.InitRand();
    m_nNodes = nNodes;
    m_nEdges = 0;
    m_fEpsilon = fEpsilon;
    m_pCorrMatrix = NULL;
    m_fAssortCoeff = -2;
    m_vNodeDegree.reserve(nNodes);
    for(int i = 1; i <= nInitNetworkSize; i++)
    {
        for(int j = 1; j <= nInitNetworkSize; j++)
        {
            if(i != j)
            {
                m_mmEdgeList.insert(pair<int ,int>(i,j));
                m_nEdges ++;
            }
        }
        m_vNodeDegree.push_back(nInitNetworkSize-1);
    }

    for(int i = nInitNetworkSize; i < nNodes; i++)
    {
        AddNewNode(i+1, nDegreeArrivingNodes);
    }

    WriteEdgeListtoFile("BA_Edge_List.txt");
    WriteNodeDegreetoFile("BA_Node_Degree.txt");

    m_nMaxDeg = FindMaxDeg();
    GenerateCorrelationMatrix();
    m_fAssortCoeff = AssortativityCoeff();

    GenerateDegDist();
}

BANetwork::~BANetwork()
{
    for(int i = 0; i <= m_nMaxDeg; i++)
    {
        delete[] m_pCorrMatrix[i];
    }

    delete[] m_pCorrMatrix;
}

void BANetwork::AddNewNode(int nNodeID, int nDegreeArrivingNodes)
{
    map<int, float> mProbLine;
    float ntemp = 0;

    for (int i = 0; i < nNodeID-1; i++)
    {
        float deg_prob = (m_vNodeDegree[i]*1.0f + m_fEpsilon)/(m_nEdges + (m_fEpsilon*(nNodeID-1)));
        mProbLine.insert(pair<int, float>(i+1, ntemp+deg_prob));
        ntemp = ntemp + deg_prob;
    }

    #ifdef DBGPRNT
    cout << "Prob line generated. ntemp = "<< ntemp << endl;
    #endif

    for (int i = 1; i <= nDegreeArrivingNodes; i++)
    {
        float prob = m_basfunc.RandNum()%100;
        prob = prob/100;

        for(map<int, float>::iterator it = mProbLine.begin(); it != mProbLine.end(); it++)
        {
            if ((*it).second > prob)
            {
                int nNodetoConnect = (*it).first;

                if(CanEdgeExist(nNodeID, nNodetoConnect))
                {
                    m_mmEdgeList.insert(pair<int ,int>(nNodetoConnect, nNodeID));
                    m_mmEdgeList.insert(pair<int ,int>(nNodeID, nNodetoConnect));
                    m_nEdges = m_nEdges + 2;

                    m_vNodeDegree[nNodetoConnect - 1] ++;
                }
                else
                {
                    i --;
                }
                break;
            }
        }
    }
    m_vNodeDegree.push_back(nDegreeArrivingNodes);
}

bool BANetwork::CanEdgeExist(int nNode1, int nNode2)
{
    for(map<int, int>::iterator it = m_mmEdgeList.begin(); it != m_mmEdgeList.end(); it ++)
    {
        if(((*it).first == nNode1 && (*it).second == nNode2) || ((*it).first == nNode2 && (*it).second == nNode2))
        {
            return false;
        }
    }

    return true;
}

int BANetwork::FindMaxDeg()
{
    int max = -1;

    for(int i = 0; i < m_vNodeDegree.size(); i++)
    {
        if(max < m_vNodeDegree[i])
        {
            max = m_vNodeDegree[i];
        }
    }

    return max;
}

void BANetwork::GenerateCorrelationMatrix()
{
    if(m_pCorrMatrix ==  NULL)
    {
        m_pCorrMatrix = new float*[m_nMaxDeg + 1];

        for(int i=0; i <= m_nMaxDeg; i++)
        {
            m_pCorrMatrix[i] = new float[m_nMaxDeg + 1];
        }

        for(int i = 0; i <= m_nMaxDeg ; i++)
        {
            for(int j = 0; j <= m_nMaxDeg ; j++)
            {
                m_pCorrMatrix[i][j] = 0;
            }
        }

        for(map<int, int>::iterator it = m_mmEdgeList.begin(); it != m_mmEdgeList.end(); it ++)
        {
            int nDeg1 = m_vNodeDegree[(*it).first - 1];
            int nDeg2 = m_vNodeDegree[(*it).second - 1];

            m_pCorrMatrix[nDeg1][nDeg2] ++;
        }

        float chk = 0;

        for(int i = 0; i <= m_nMaxDeg ; i++)
        {
            for(int j = 0; j <= m_nMaxDeg ; j++)
            {
                m_pCorrMatrix[i][j] = m_pCorrMatrix[i][j]/m_nEdges;
                chk = chk + m_pCorrMatrix[i][j];
            }
        }

        cout << "chk:: " << chk << endl;
    }
    else
    {
        cout << "Correlation Matrix already exists" << endl;
    }
}

float BANetwork::AssortativityCoeff()
{
    if( m_pCorrMatrix == NULL)
    {
        cout << "The correlation matrix doesn't exist" << endl;
        return -2;
    }
    else
    {
        float trm1 = 0.0f, trm2 = 0.0f, trm3 = 0.0f, r = 0.0f;

        for(int i = 0; i <= m_nMaxDeg; i++)
        {
            for(int j = 0; j <= m_nMaxDeg; j++)
            {
                float fac = m_pCorrMatrix[i][j];
                trm1 = trm1 + fac*i*j;
                trm2 = trm2 + fac*((i+j)/2.0);
                trm3 = trm3 + fac*(((i*i)+(j*j))/2.0);
            }
        }

        trm2 = trm2*trm2;

        r = (trm1 - trm2)/(trm3 - trm2);
        return r;
    }
}

void BANetwork::GenerateDegDist()
{
    m_vDegDist.erase(m_vDegDist.begin(), m_vDegDist.end());
    for(int i = 0; i <= m_nMaxDeg; i++)
    {
        m_vDegDist.push_back(0);
    }

    for(int i = 0; i < m_vNodeDegree.size(); i++)
    {
        if(m_vNodeDegree[i] >= 0)
        {
            m_vDegDist[m_vNodeDegree[i]] ++;
        }
    }

    for(int i = 0; i < m_vDegDist.size(); i++)
    {
        m_vDegDist[i] = (m_vDegDist[i]*1.0f)/m_nNodes;
    }
}
void BANetwork::ApplyAttack(int choice, float del_fr)
{
    int nNumofNodesDel = del_fr * m_nNodes, count = 0;

    vector<int> vRemainingNodeID;
    vector<int> vDelNodeID;

    for(int i = 0; i < m_nNodes; i++)
    {
        vRemainingNodeID.push_back(0);
    }

    if( choice == RANDOM_FAILURE)
    {
        while(count < nNumofNodesDel)
        {
            int sel = m_basfunc.RandNum() % m_nNodes;
            if( vRemainingNodeID[sel] != -1)
            {
                count ++;
                vRemainingNodeID[sel] = -1;
                vDelNodeID.push_back(sel+1);
            }
        }
    }

    int nEdgeLiSize = m_mmEdgeList.size();
    int **pEdge = new int*[nEdgeLiSize];

    for(int i=0; i < nEdgeLiSize; i++)
    {
        pEdge[i] = new int[2];
    }

    int edgecount = 0;
    for(map<int, int>::iterator it = m_mmEdgeList.begin(); it != m_mmEdgeList.end(); it ++)
    {
        pEdge[edgecount][0] = (*it).first;
        pEdge[edgecount][1] = (*it).second;
        edgecount ++;
    }

    for(int i = 0; i < edgecount; i++)
    {
        int nNode1 = pEdge[i][0];
        int nNode2 = pEdge[i][1];

        if(vRemainingNodeID[nNode1-1] == -1 || vRemainingNodeID[nNode2-1] == -1)
        {
            //m_mmEdgeList.erase(it);
            pEdge[i][0] = pEdge[i][1] = -1;
            m_vNodeDegree[nNode1-1] --;
            m_nEdges --;
        }
    }
    /*for(map<int, int>::iterator it = m_mmEdgeList.begin(); it != m_mmEdgeList.end(); it ++)
    {
        int nNode1 = (*it).first;
        int nNode2 = (*it).second;

        if(vRemainingNodeID[nNode1-1] == -1 || vRemainingNodeID[nNode2-1] == -1)
        {
            m_mmEdgeList.erase(it);
            m_vNodeDegree[nNode1-1] --;
            m_nEdges --;
        }
    }*/

    for(int i = 0; i < vDelNodeID.size(); i++)
    {
        if( m_vNodeDegree[vDelNodeID[i] - 1] != 0)
        {
            cout << "Something is wrong here" << endl;
        }
        else
        {
            m_vNodeDegree[vDelNodeID[i] - 1] = -1;
        }
    }

    m_mmEdgeList.erase(m_mmEdgeList.begin(), m_mmEdgeList.end());
    for(int i = 0; i < edgecount; i++)
    {
        if( pEdge[i][0] != -1)
        {
            m_mmEdgeList.insert(pair<int ,int>(pEdge[i][0], pEdge[i][1]));
        }
    }

    for(int i = 0; i < edgecount; i++)
    {
        delete[] pEdge[i];
    }
    delete[] pEdge;
    pEdge = NULL;

    m_nNodes = m_nNodes - nNumofNodesDel;

    for(int i = 0; i <= m_nMaxDeg; i++)
    {
        delete[] m_pCorrMatrix[i];
    }
    delete[] m_pCorrMatrix;
    m_pCorrMatrix = NULL;

    WriteEdgeListtoFile("BA_Edge_List_After_Attack.txt");
    WriteNodeDegreetoFile("BA_Node_Degree_After_Attack.txt");

    m_nMaxDeg = FindMaxDeg();
    GenerateCorrelationMatrix();
    m_fAssortCoeff = AssortativityCoeff();

    GenerateDegDist();
}
void BANetwork::ApplyAttack_deterministically(int choice,float del_fr)
{

    int nNumofNodesDel = del_fr * m_nNodes, count = 0;

    vector<int> vRemainingNodeID;
    vector<int> vDelNodeID;

    for(int i = 0; i < m_nNodes; i++)
    {
        vRemainingNodeID.push_back(0);
    }
    bool flag=false;
    vector<int> node_degree1;
    if( choice == 2)
    {

        int z=0;

        while(count < nNumofNodesDel)
        {
           int sel;
            vector<int> node_degree;
            for(int i=0;i<m_vNodeDegree.size();i++){
                if(flag==false)
                  flag=true;
                node_degree.push_back(m_vNodeDegree[i]);
                node_degree1.push_back(m_vNodeDegree[i]);
            }
           sort(node_degree.begin(),node_degree.end());
            int x=node_degree[node_degree.size()-1];
            vector<int>::iterator it=m_vNodeDegree.begin();
            for(int i=0;i<m_vNodeDegree.size();i++,it++)
            {
                if(x==m_vNodeDegree[i])
                    {
                        sel=i;
                        m_vNodeDegree.erase(it);
                        break;
                    }

            }
            node_degree.erase(node_degree.end()-1);
            if( vRemainingNodeID[sel] != -1)
            {
                count++;
                vRemainingNodeID[sel] = -1;
                vDelNodeID.push_back(sel+1);

        }

    }
  }
    if(flag==true)
    m_vNodeDegree.swap(node_degree1);

    int nEdgeLiSize = m_mmEdgeList.size();
    int **pEdge = new int*[nEdgeLiSize];

    for(int i=0; i < nEdgeLiSize; i++)
    {
        pEdge[i] = new int[2];
    }

    int edgecount = 0;
    for(map<int, int>::iterator it = m_mmEdgeList.begin(); it != m_mmEdgeList.end(); it ++)
    {
        pEdge[edgecount][0] = (*it).first;
        pEdge[edgecount][1] = (*it).second;
        edgecount ++;
    }

    for(int i = 0; i < edgecount; i++)
    {
        int nNode1 = pEdge[i][0];
        int nNode2 = pEdge[i][1];

        if(vRemainingNodeID[nNode1-1] == -1 || vRemainingNodeID[nNode2-1] == -1)
        {
            //m_mmEdgeList.erase(it);
            pEdge[i][0] = pEdge[i][1] = -1;
            m_vNodeDegree[nNode1-1] --;
            m_nEdges --;
        }
    }
    /*for(map<int, int>::iterator it = m_mmEdgeList.begin(); it != m_mmEdgeList.end(); it ++)
    {
        int nNode1 = (*it).first;
        int nNode2 = (*it).second;

        if(vRemainingNodeID[nNode1-1] == -1 || vRemainingNodeID[nNode2-1] == -1)
        {
            m_mmEdgeList.erase(it);
            m_vNodeDegree[nNode1-1] --;
            m_nEdges --;
        }
    }*/

    for(int i = 0; i < vDelNodeID.size(); i++)
    {
        if( m_vNodeDegree[vDelNodeID[i] - 1] != 0)
        {
            cout << "Something is wrong here" << endl;
        }
        else
        {
            m_vNodeDegree[vDelNodeID[i] - 1] = -1;
        }
    }
    m_mmEdgeList.erase(m_mmEdgeList.begin(), m_mmEdgeList.end());
    for(int i = 0; i < edgecount; i++)
    {
        if( pEdge[i][0] != -1)
        {
            m_mmEdgeList.insert(pair<int ,int>(pEdge[i][0], pEdge[i][1]));
        }
    }

    for(int i = 0; i < edgecount; i++)
    {
        delete[] pEdge[i];
    }
    delete[] pEdge;
    pEdge = NULL;

    m_nNodes = m_nNodes - nNumofNodesDel;
    for(int i = 0; i <= m_nMaxDeg; i++)
    {
        delete[] m_pCorrMatrix[i];
    }
    delete[] m_pCorrMatrix;
    m_pCorrMatrix = NULL;
    WriteEdgeListtoFile("BA_Edge_List_After_Attack.txt");
    WriteNodeDegreetoFile("BA_Node_Degree_After_Attack.txt");

    m_nMaxDeg = FindMaxDeg();
    GenerateCorrelationMatrix();
    m_fAssortCoeff = AssortativityCoeff();

    GenerateDegDist();

}
void BANetwork::WriteEdgeListtoFile(char *fname)
{
    fstream fwr(fname, ios::out);

    for(map<int, int>::iterator it = m_mmEdgeList.begin(); it != m_mmEdgeList.end(); it ++)
    {
        fwr << (*it).first << " " << (*it).second <<endl;
    }
    fwr.close();
}
void BANetwork::WriteNodeDegreetoFile(char *fname)
{
    fstream fwr(fname, ios::out);

    for(int i = 0; i < m_vNodeDegree.size(); i++)
    {
        fwr << (i+1) << " " << m_vNodeDegree[i] << endl;
    }
    fwr.close();
}

void BANetwork::WriteDegDisttoFile(char *fname)
{
    fstream fwr(fname, ios::out);

    for(int i = 0; i < m_vDegDist.size(); i++)
    {
        fwr << i << " " << m_vDegDist[i] << endl;
    }
    fwr.close();
}
float BANetwork::GetAssortCoeff()
{
    return m_fAssortCoeff;
}
vector<float>& BANetwork::GetDegDist()
{
    return m_vDegDist;
}
void BANetwork::degree_distribution(char *fname)
{
  int a[10000];
  fstream f(fname,ios::out);
  fill_n(a,10000,0);
  for(int i=0;i<m_vNodeDegree.size();i++)
  {
        a[m_vNodeDegree[i]]++;
  }
  cout<<"____________________________________________________\n";
  for(int i=0;i<m_vNodeDegree.size();i++)
  {
        if(a[i]>0)
        cout<<"count for degree :"<<i<<" "<<a[i]<<endl;
  }

  vector<int> pk;
  int i1=1;
  int i2=2;
  f<<i1<<" "<<i2<<endl;
 for(int i=0;i<10000;i++)
  {
     if(a[i]>0){
        int k=a[i];
        int pp=m_vNodeDegree.size();
        float c=(float)k/(float)pp;
        f<<i<<" "<<c<<endl;

     }

   }
  f.close();
}
/*
int BANetwork::factorial(int x)
{
    int fac=1;
    for(int i=x;i>0;i--)
        fac=fac*i;
    return fac;

}
int  BANetwork::combination(int m,int n)
{
    return factorial(m)/(factorial(m-n)*factorial(n));

}
void BANetwork::clustering_coefficient()
{
    int sum=0;
    for(int i=0;i<m_vNodeDegree.size();i++)
    {
        vector<int> neighbour;
        for(int j=0;j<m_vNodeDegree.size();j++)
        {
            if(i!=j)
            {
                bool flag = CanEdgeExist(i,j);
                if(flag==false)
                    neighbour.push_back(j);
            }
        }
        int n=0;
        for(vector<int>::iterator it=neighbour.begin();it!=neighbour.end();it++)
        {

                for(vector<int>::iterator itt=neighbour.begin();itt!=neighbour.end();itt++)
                   {
                       if(*it!=*itt)
                          {
                              bool flag = CanEdgeExist(*it,*itt);
                              if(flag==true)
                                n++;
                          }
                   }
        }
        n=n/2;
        int p=m_vNodeDegree[i];
        int denominator=combination(p,2);
        sum+=n/denominator;
    }
    int clustering = sum/m_vNodeDegree.size();
    cout<<"The clustring coefficient of network is : "<<clustering<<endl;


}
/*
void BANetwork::findlength(int source)
 {
     vector<int> dist;

     for(map<int,int>::iterator it=m_mmEdgeList.begin();it!=m_mmEdgeList.end();it++)
     {
         set<int> neighbour;
         if((*it).first==source||(*it).second==source)



     }


 }
void BANetwork::diameter()
{
    int max=INT_MIN;
    for(int i=0;i<m_vNodeDegree.size();i++)
    {
         findlength(i);
    }



}
*/
int main()
{

   int nTotalNodes = 0, nInitNetworkSize = 0, nDegreeArrival = 0, choice = INT_MIN;
    float fEpsilon = 0.0f, del_fr = 0.0f;

    cout << "Enter total no. of nodes:: ";
    cin >> nTotalNodes;

    cout << "Enter initial complete network size:: ";
    cin >> nInitNetworkSize;

    cout << "Enter the degree of the arriving nodes:: ";
    cin >> nDegreeArrival;

    cout << "Enter epsilon:: ";
    cin >> fEpsilon;

    BANetwork BANet(nTotalNodes, nInitNetworkSize, nDegreeArrival, fEpsilon);

    cout << "Network generation done" << endl;
    cout<<"degree distribution\n";
    BANet.degree_distribution("degree_distribution_before_attack.txt");

    cout << "Assortativity coeff of the network:: " << BANet.GetAssortCoeff() << endl;
    BANet.WriteDegDisttoFile("BA_Deg_Dist.txt");

    cout << "\nRandom failure ------------> 1" << endl;
    cout << "Deterministic attack ------> 2" << endl;
    cout<<"Exit---------------------------- -1"<<endl;
    cout << "Enter your choice:: " ;
    cin >> choice;
   char  fname[100];
   switch(choice)
    {
        case 1:
        cout <<"Enter fraction of nodes to be removed randomly:: ";
        cin >> del_fr;
        BANet.ApplyAttack(choice, del_fr);
        cout<<"Enter the name of file : ";
        cin>>fname;
        BANet.degree_distribution(fname);
        break;

        case 2:
        cout <<"Enter fraction of nodes to be removed deterministically:: ";
        cin >> del_fr;
        BANet.ApplyAttack_deterministically(choice, del_fr);
        break;

        default:
        cout << "Invalid choice" << endl;
    }


      Graph g;
      multimap<int,int>::iterator it;
      for(it=BANet.m_mmEdgeList.begin();it!=BANet.m_mmEdgeList.end();it++){
        int x=(*it).first;
        int y=(*it).second;
        g.addEdge(x,y);
      }
         cout << "Following are strongly connected components in "
            "given graph \n";
           g.printSCCs(BANet,del_fr);
      //  cout<<maximum<<endl;




        cout << "Attack simulation done" << endl;
        cout << "Assortativity coeff of the network:: " << BANet.GetAssortCoeff() << endl;

        char file[]="BA_Deg_Dist_After_Attack.txt";

        BANet.WriteDegDisttoFile(file);



}
