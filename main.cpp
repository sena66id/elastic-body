#include<iostream>
#include<fstream>
#include<vector>
#include<string>
#include<sstream>
#include<cmath>
#include "Eigen/Core"
#include "Eigen/LU"
#define E 205000.0//N/mm*2
#define nu 0.3
#include "Eigen/SparseCore"
#include "Eigen/SparseLU"

double calcDeterminant_3x3(const double (&a)[3][3]);
void calcInverseMatrix_3x3(double (&inv_a)[3][3],const double (&a)[3][3]);



using namespace Eigen;
using namespace std;

typedef Triplet<double> T;


void export_vtu(const std::string &file, vector<vector<double>> node, vector<vector<int>> element, vector<double> u)
{
    FILE *fp;
  if ((fp = fopen(file.c_str(), "w")) == NULL)
  {
    cout << file << " open error" << endl;
    exit(1);
  }

  fprintf(fp, "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt32\">\n");
  fprintf(fp, "<UnstructuredGrid>\n");
  fprintf(fp, "<Piece NumberOfPoints= \"%d\" NumberOfCells= \"%d\" >\n", node.size(), element.size());
  fprintf(fp, "<Points>\n");
  int offset = 0;
  fprintf(fp, "<DataArray type=\"Float64\" Name=\"Position\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%d\"/>\n",offset);
  offset += sizeof(int) + sizeof(double) * node.size() * 3;
  fprintf(fp, "</Points>\n");

  fprintf(fp, "<Cells>\n");
  fprintf(fp, "<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n");
  for (int i = 0; i < element.size(); i++){
    for (int j = 0; j < element[i].size(); j++) fprintf(fp, "%d ", element[i][j]);

    
    fprintf(fp, "\n");
  }
  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n");
  int num = 0;
  for (int i = 0; i < element.size(); i++)
  {
    num += element[i].size();
    fprintf(fp, "%d\n", num);
  }
  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
  for (int i = 0; i < element.size(); i++) fprintf(fp, "%d\n", 10);
    
  fprintf(fp, "</DataArray>\n");
  fprintf(fp, "</Cells>\n");

  fprintf(fp, "<PointData>\n");

  fprintf(fp, "<DataArray type=\"Float64\" Name=\"displacement[m]\" NumberOfComponents=\"3\" format=\"appended\" offset=\"%d\"/>\n",offset);
  offset += sizeof(int) + sizeof(double) * node.size()*3;

  fprintf(fp, "</PointData>\n");

  fprintf(fp, "<CellData>\n");
  fprintf(fp, "</CellData>\n");
  fprintf(fp, "</Piece>\n");
  fprintf(fp, "</UnstructuredGrid>\n");
  fprintf(fp, "<AppendedData encoding=\"raw\">\n");
  fprintf(fp, "_");
  fclose(fp);

  fstream ofs;
  ofs.open(file.c_str(), ios::out | ios::app | ios_base::binary);
  double *data_d = new double[node.size()*3];
  num = 0;
  int size=0;
  for (int ic = 0; ic < node.size(); ic++){
    data_d[num] = node[ic][0];
    num++;
    data_d[num] = node[ic][1];
    num++;
    data_d[num] = node[ic][2];
    num++;
  }
  size=sizeof(double)*node.size()*3;
  ofs.write((char *)&size, sizeof(size));
  ofs.write((char *)data_d, size);

  num=0;
  for (int ic = 0; ic < node.size(); ic++){
      data_d[num]   = u[ic];
      num++;
      data_d[num]   = u[node.size()+ic];
      num++;
      data_d[num]   = u[node.size()*2+ic];
      num++;
  }
  size=sizeof(double)*node.size()*3;
  ofs.write((char *)&size, sizeof(size));
  ofs.write((char *)data_d, size);

  delete data_d;

  ofs.close();

  if ((fp = fopen(file.c_str(), "a")) == NULL)
  {
    cout << file << " open error" << endl;
    exit(1);
  }
  fprintf(fp, "\n</AppendedData>\n");
  fprintf(fp, "</VTKFile>\n");
  fclose(fp);
}

int main()
{

    string str;
    ifstream ifs("node.dat");
    //vector<double> t(2) は double t[2]と同じ
    //vector<vector<double>> t(2, vector<double>(2))は double t[2][2]と同じ
    
    vector<vector<double>> x;
    while(getline(ifs,str)){
        istringstream ss(str);
        string tmp;
        vector<double> tmp_x;
        for(int j=0; j<3; j++){
            getline(ss, tmp, ' ');
            tmp_x.push_back(stod(tmp));
        }
        x.push_back(tmp_x);
    }
    ifs.close();

    ifs.open("element.dat");
    vector<vector<int>> element;
    while(getline(ifs,str)){
        istringstream ss(str);
        string tmp;
        vector<int> tmp_element;
        for(int j=0; j<4; j++){
            getline(ss, tmp, ' ');
            tmp_element.push_back(stoi(tmp));
        }
        element.push_back(tmp_element);
    }
    ifs.close();

    vector<double> C1(x.size(),0.0);

    ofstream ofs1("boundary_left.dat");
    for(int i=0; i<x.size(); i++){
        if(fabs(x[i][2]-0.0)<0.000001){
            ofs1 << i << endl;
            C1[i] = 1.0;
        }
    }
    ofs1.close();

    vector<double> C2(x.size(),0.0);
    ofstream ofs2("boundary_right.dat");
    for(int i=0; i<x.size(); i++){
        if(fabs(x[i][2]-0.0)==0.1){
            ofs2 << i << endl;
            C2[i] = 1.0;
        }
    }
    ofs2.close();

    //export_vtu("test.vtu", x, element, C1);

     ifs.open("boundary_left.dat");
     vector<int> boundary_left;
     while(getline(ifs,str)){
        istringstream ss(str);
        string tmp;
        getline(ss, tmp, ' ');
          boundary_left.push_back(stoi(tmp));
      }
      ifs.close();

      ifs.open("boundary_right.dat");
      vector<int> boundary_right;
      while(getline(ifs,str)){
          istringstream ss(str);
          string tmp;
          getline(ss, tmp, ' ');
          boundary_right.push_back(stoi(tmp));
      }
      ifs.close();



    MatrixXd K(x.size()*3, x.size()*3);
    MatrixXd K2(x.size()*3, x.size()*3);
    VectorXd U(x.size()*3);
    VectorXd F(x.size()*3);
    VectorXd R(x.size()*3);
    U = VectorXd::Zero(x.size()*3);
    F = VectorXd::Zero(x.size()*3);
    R = VectorXd::Zero(x.size()*3);
    K = MatrixXd::Zero(x.size()*3,x.size()*3);

    
    double l=E*nu/((1.0+nu)*(1.0-2.0*nu));
    double mew=E/(2.0*(1.0+nu));
    vector<vector<double>> D{
      {l+2.0*mew,l,l,0,0,0},
      {l,l+2.0*mew,l,0,0,0},
      {l,l,l+2.0*mew,0,0,0},
      {0,0,0,mew,0,0},
      {0,0,0,0,mew,0},
      {0,0,0,0,0,mew}
    };
    
    double dNdr[3][4]; 
    dNdr[0][0] = -1.0; dNdr[0][1] = 1.0; dNdr[0][2] = 0.0; dNdr[0][3] = 0.0;
    dNdr[1][0] = -1.0; dNdr[1][1] = 0.0; dNdr[1][2] = 1.0; dNdr[1][3] = 0.0;
    dNdr[2][0] = -1.0; dNdr[2][1] = 0.0;  dNdr[2][2] = 0.0; dNdr[2][3] = 1.0;
    

    double dxdr[3][3];
    //3*x.size(),vector<double>(x.size()*3,0.0));

    for(int ic=0; ic<element.size(); ic++){
      //cout << ic << endl;
      for(int i=0; i<3; i++){
        for(int j=0; j<3; j++){
          dxdr[i][j]=0.0;
          for(int k=0; k<4; k++){
            dxdr[i][j] += dNdr[i][k]*x[element[ic][k]][j];
          }
        }
      }

      
      
      double dxdr2[3][3];

       for(int i=0; i<3; i++){
        for(int j=0; j<3; j++){
            dxdr2[i][j]=dxdr[j][i];
        }
       }

      


      double inv_dxdr2[3][3];
      double det;


      det = calcDeterminant_3x3(dxdr);
      calcInverseMatrix_3x3(inv_dxdr2,dxdr);

      double dndx[3][4];
      
      
      for(int i=0;i<3;i++){
            for(int j=0;j<4;j++){
                dndx[i][j]=0;
                for(int k=0;k<3;k++){
                    dndx[i][j]+=inv_dxdr2[i][k]*dNdr[k][j];
                }
            }
     }

      double B[6][12];
      
      double Bn[6][12];
      
      for(int a=4; a<12;a++){
        B[0][a]=0;
        }
      

      for(int a=0; a<12;a++){
        if(a!=4||a!=5||a!=6||a!=7){
          B[1][a]=0;
        }
      }

      for(int a=0; a<8;a++){
        B[2][a]=0;
        }
      for(int a=0; a<4; a++){
        B[3][a]=0;
      }
      for(int a=0; a<12;a++){
        if(a==4 || a==5 || a==6||a==7){
          B[4][a]=0;
        }
      }
      for(int a=8; a<12;a++){
          B[5][a]=0;
        }
      

      B[0][0]=dndx[0][0]; B[0][1]=dndx[0][1]; B[0][2]=dndx[0][2]; B[0][3]=dndx[0][3];
      B[1][4]=dndx[1][0]; B[1][5]=dndx[1][1]; B[1][6]=dndx[1][2]; B[1][7]=dndx[1][3];
      B[2][8]=dndx[2][0]; B[2][9]=dndx[2][1]; B[2][10]=dndx[2][2]; B[2][11]=dndx[2][3];
      B[3][4]=dndx[2][0]; B[3][5]=dndx[2][1]; B[3][6]=dndx[2][2]; B[3][7]=dndx[2][3]; B[3][8]=dndx[1][0]; B[3][9]=dndx[1][1]; B[3][10]=dndx[1][2]; B[3][11]=dndx[1][3];
      B[4][0]=dndx[2][0]; B[4][1]=dndx[2][1]; B[4][2]=dndx[2][2]; B[4][3]=dndx[2][3]; B[4][8]=dndx[0][0]; B[4][9]=dndx[0][1]; B[4][10]=dndx[0][2]; B[4][11]=dndx[0][3];
      B[5][0]=dndx[1][0]; B[5][1]=dndx[1][1]; B[5][2]=dndx[1][2]; B[5][3]=dndx[1][3]; B[5][4]=dndx[0][0]; B[5][5]=dndx[0][1]; B[5][6]=dndx[0][2]; B[5][7]=dndx[0][3];
           
      double trans_B[12][6];
      double r1[12][6];
      double Ke[12][12];

      for(int g=0;g<12;g++){
        for(int h=0;h<6;h++){
          trans_B[g][h]=B[h][g];
          //cout<<trans_B[g][h]<<endl;
        }
      }
          
      for(int i=0;i<12;i++){
        for(int j=0; j<6; j++){
          r1[i][j]=0;
          for(int k=0; k<6; k++){
            r1[i][j]+=trans_B[i][k]*D[k][j];
          }
        }
      }
      for(int i=0;i<12;i++){
        for(int j=0; j<12; j++){
          Ke[i][j]=0;
          for(int k=0; k<6; k++){
            Ke[i][j]+=r1[i][k]*B[k][j]*calcDeterminant_3x3(dxdr2)*0.167;
          }
        }
      }
      

     for(int i=0;i<4;i++){
        for(int j=0;j<4;j++){
            K(element[ic][i],element[ic][j]) += Ke[i][j];
        }

     }
     for(int i=0;i<4;i++){
        for(int j=4;j<8;j++){
            K(element[ic][i],element[ic][j-4]+x.size()) += Ke[i][j];
        }
    }
     for(int i=0;i<4;i++){
        for(int j=8;j<12;j++){
            K(element[ic][i],element[ic][j-8]+2*x.size()) += Ke[i][j];
        }
    }
     for(int i=4;i<8;i++){
        for(int j=0;j<4;j++){
            K(element[ic][i-4]+x.size(),element[ic][j])+= Ke[i][j];
        }

     }
     for(int i=4;i<8;i++){
        for(int j=4;j<8;j++){
            K(element[ic][i-4]+x.size(),element[ic][j-4]+x.size()) += Ke[i][j];
        }
    }
     for(int i=4;i<8;i++){
        for(int j=8;j<12;j++){
            K(element[ic][i-4]+x.size(),element[ic][j-8]+2*x.size()) += Ke[i][j];
        }
    }
    for(int i=8;i<12;i++){
        for(int j=0;j<4;j++){
            K(element[ic][i-8]+2*x.size(),element[ic][j]) += Ke[i][j];
        }

     }
     for(int i=8;i<12;i++){
        for(int j=4;j<8;j++){
            K(element[ic][i-8]+2*x.size(),element[ic][j-4]+x.size()) += Ke[i][j];
        }
    }
     for(int i=8;i<12;i++){
        for(int j=8;j<12;j++){
            K(element[ic][i-8]+2*x.size(),element[ic][j-8]+2*x.size()) += Ke[i][j];
        }
    }
    }

    for(int i=0;i<x.size()*3;i++){
      for(int j=0; j<x.size()*3;j++){
        K2(i,j)=K(i,j);
        
      }
    }

    

    for(int i=0; i<boundary_left.size(); i++){
      for(int j=0; j<x.size()*3; j++){
        K(boundary_left[i],j) = 0;
        K(boundary_left[i]+x.size(),j) = 0;
        K(boundary_left[i]+2*x.size(),j) = 0;
      }
    }

    for(int i=0; i<boundary_right.size(); i++){
      for(int j=0; j<x.size()*3; j++){
      //K(boundary_right[i],j) = 0;
      //K(boundary_right[i]+x.size(),j) = 0;
      K(boundary_right[i]+2*x.size(),j) = 0;
      }
    }

    for(int i=0; i<boundary_left.size(); i++){
        K(boundary_left[i],boundary_left[i]) = 1.0;
        K(boundary_left[i]+x.size(),boundary_left[i]+x.size()) = 1.0;
        K(boundary_left[i]+2*x.size(),boundary_left[i]+2*x.size()) = 1.0;
      }

      for(int i=0; i<boundary_right.size(); i++){
        //K(boundary_right[i],boundary_right[i]) = 1.0;
       //K(boundary_right[i]+x.size(),boundary_right[i]+x.size()) = 1.0;
        K(boundary_right[i]+2*x.size(),boundary_right[i]+2*x.size()) = 1.0;
      }

    for(int i=0; i<boundary_right.size(); i++){
        R(boundary_right[i]+2*x.size()) = 0.00001;
      }


      

      SparseMatrix<double> K_sparse(x.size()*3, x.size()*3);
      vector<T> tripletList;
      for (int i=0;i<x.size()*3 ;i++){
        for(int j=0;j<x.size()*3 ;j++){
          if(K(i,j)!=0){
            tripletList.push_back(T(i,j,K(i,j)));
          }

        }
      }
      K_sparse.setFromTriplets(tripletList.begin(),tripletList.end());

     SparseLU<SparseMatrix<double>, COLAMDOrdering<int>> solver;
    solver.compute(K_sparse);
    if(solver.info()!=Success) {
    // decomposition failed
    cout << "decomposition failed" << endl;
    return -1;
    }
    U = solver.solve(R);
    if(solver.info()!=Success) {
    // solving failed
    std::cout << "solving failed" << endl;
    return -1;
    }
    
    for(int i=0; i<x.size(); i++){
        x[i][0] += U[i];
        x[i][1] += U[i+x.size()];
        x[i][2] += U[i+x.size()*2];

        
    }

  

    for(int i=0;i<x.size()*3;i++){
      for(int k=0; k<x.size()*3;k++){
     
          F(i)+=K2(i,k)*U(k);
        }
      }
    ofstream outputfile("F.dat");
    outputfile << F;
    outputfile.close();

    vector<double> F1(x.size()*3);


    for(int i=0; i<x.size()*3; i++){
        F1[i]=F[i];
    }

    double sum1=0;
    for(int i=0; i<boundary_left.size();i++){
      sum1+=F1[2*x.size()+i];
    }

    //cout<<sum1<<endl;

     double sum2=0;
    for(int i=boundary_left.size(); i<boundary_right.size()+boundary_left.size();i++){
      sum2+=F1[2*x.size()+i];
    }

   

    double keisan,kaiseki;

    keisan=sum2/(0.01*0.01*3.14);

   cout<<keisan<<endl;





    

    //export_vtu("result.vtu", x, element, F1);

}

double calcDeterminant_3x3( const double (&a)[3][3])
{
  double det  = a[0][0] * a[1][1] * a[2][2] + a[1][0] * a[2][1] * a[0][2] + a[2][0] * a[0][1] * a[1][2]
              - a[2][0] * a[1][1] * a[0][2] - a[1][0] * a[0][1] * a[2][2] - a[0][0] * a[2][1] * a[1][2];
  return det;
}

void calcInverseMatrix_3x3(double (&inv_a)[3][3],const double (&a)[3][3])
{
  double det;

  det = calcDeterminant_3x3(a);

  inv_a[0][0] = a[1][1]*a[2][2] - a[1][2]*a[2][1];
  inv_a[0][1] = a[0][2]*a[2][1] - a[0][1]*a[2][2];
  inv_a[0][2] = a[0][1]*a[1][2] - a[0][2]*a[1][1];
  inv_a[1][0] = a[1][2]*a[2][0] - a[1][0]*a[2][2];
  inv_a[1][1] = a[0][0]*a[2][2] - a[0][2]*a[2][0];
  inv_a[1][2] = a[0][2]*a[1][0] - a[0][0]*a[1][2];
  inv_a[2][0] = a[1][0]*a[2][1] - a[1][1]*a[2][0];
  inv_a[2][1] = a[0][1]*a[2][0] - a[0][0]*a[2][1];
  inv_a[2][2] = a[0][0]*a[1][1] - a[0][1]*a[1][0];

  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++) inv_a[i][j] = inv_a[i][j] / det;
}
}
   // vector<vector<double>> K(


