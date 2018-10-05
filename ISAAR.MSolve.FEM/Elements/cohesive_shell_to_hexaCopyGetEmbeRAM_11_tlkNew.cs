using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.FEM.Interfaces;//using ISAAR.MSolve.PreProcessor.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;//using ISAAR.MSolve.Matrices.Interfaces;
using System.Runtime.InteropServices;
using ISAAR.MSolve.Numerical.LinearAlgebra;//using ISAAR.MSolve.Matrices;
using ISAAR.MSolve.FEM.Elements.SupportiveClasses;//using ISAAR.MSolve.PreProcessor.Elements.SupportiveClasses;
using ISAAR.MSolve.FEM.Embedding;//using ISAAR.MSolve.PreProcessor.Embedding;
// compa
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Materials.Interfaces;
using ISAAR.MSolve.FEM;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Interfaces;
using IEmbeddedElement = ISAAR.MSolve.FEM.Interfaces.IEmbeddedElement;
using ISAAR.MSolve.FEM.Interpolation;
using ISAAR.MSolve.Discretization.Integration.Quadratures;


//TODO: current nodal coordinates should be managed by the analyzer, instead of each element calculting and storing them independently. The same applies for direction vectors of shells. 
//TODO: direction vectors creation and update could be handled by a dedicated class that will be composed into this element. Which element would update them then?
//TODO: perhaps separate the cases fot when the shell is over or under the cohesive element ibto 2 subclasses
namespace ISAAR.MSolve.FEM.Elements
{
    public class cohesive_shell_to_hexaCopyGetEmbeRAM_11_tlkNew : IStructuralFiniteElement, IEmbeddedElement
    {
        //metavlhtes opws sto hexa8
        protected readonly static DOFType[] nodalDOFTypes = new DOFType[] { DOFType.X, DOFType.Y, DOFType.Z };
        protected readonly static DOFType[] nodalDOFTypes2 = new DOFType[] { DOFType.X, DOFType.Y, DOFType.Z, DOFType.RotX, DOFType.RotY };
        protected readonly static DOFType[][] dofTypes = new DOFType[][] { nodalDOFTypes2, nodalDOFTypes2, nodalDOFTypes2,
            nodalDOFTypes2, nodalDOFTypes2, nodalDOFTypes2, nodalDOFTypes2, nodalDOFTypes2,nodalDOFTypes, nodalDOFTypes, nodalDOFTypes,
            nodalDOFTypes, nodalDOFTypes, nodalDOFTypes, nodalDOFTypes, nodalDOFTypes };
        protected readonly ICohesiveZoneMaterial3D[] materialsAtGaussPoints;
        private readonly InterpolationShell8CohesiveNew interpolation = InterpolationShell8CohesiveNew.UniqueInstance;
        protected IElementDOFEnumerator dofEnumerator = new GenericDOFEnumerator();
        // ews edw

        public int gp_d1_coh { get; set; } // den prepei na einai static--> shmainei idio gia ola taantikeimena afthw ths klashs
        public int gp_d2_coh { get; set; }
        private int nGaussPoints;
        public double[][] oVn_i { get; set; } // apo to shell8disp kai oti allo xreiazetai tha paroume apo shell8disp kai ena bool akomh
        public double[] tk { get; set; } //
        /// <summary>
        /// zero denotes bottom surface of the shell element
        /// </summary>
        public int ShellElementSide { get; set; } // 

        // diadikasies dhmiourgias kai ananewshs twn dianusmatwn katefthunshs apo to shell8disp kai aparaithtes metavlhtes aftwn
        private double[][] tU;   //8 arrays twn 6 stoixeiwn 
        private double[][] tUvec;//8 arrays twn 6 stoixeiwn
        //private double tV1norm;

        /// <summary>
        /// Initial coordinates of shell midsurface nodes
        /// </summary>
        private double[][] ox_i_shell_midsurface;
        //private double[][] tx_i_shell_midsurface;

        // metavlhtes vasikwn methodwn apo to cohesive 16 node

        /// <summary>
        /// Initial nodel coordinates of 16 node inner cohesive element
        /// </summary>
        private double[][] ox_i; //den einai apo afta pou orizei o xrhsths

        /// <summary>
        /// Unrolled current nodal coordinates of 16 node inner cohesive element
        /// </summary>
        private double[] x_local; // to dianusma x ths matlab sunarthshs pou einai apo t_x_global_pr

        /// <summary>
        ///  ananewsh strofwn tou stoixeiou
        /// voithitikes metavlhtes gia upologismo strofhs
        /// </summary>      
        private double[] ak_total = new double[8];
        private double[] bk_total = new double[8];
        public bool MatrixIsNotInitialized = true;


        protected cohesive_shell_to_hexaCopyGetEmbeRAM_11_tlkNew()//consztructor apo to hexa8
        {
        }

        public cohesive_shell_to_hexaCopyGetEmbeRAM_11_tlkNew(ICohesiveZoneMaterial3D material, int gp_d1c, int gp_d2c, IQuadrature2D quadratureForStiffness)
        {

            this.gp_d1_coh = gp_d1c;
            this.gp_d2_coh = gp_d2c;
            this.nGaussPoints = this.gp_d1_coh * this.gp_d2_coh;
            materialsAtGaussPoints = new ICohesiveZoneMaterial3D[nGaussPoints];
            for (int i = 0; i < nGaussPoints; i++)
                materialsAtGaussPoints[i] = material.Clone();
            this.QuadratureForStiffness = quadratureForStiffness;
        }

        public IQuadrature2D QuadratureForStiffness { get;  }
        
        private void GetInitialGeometricDataForMidsurface(IElement element) //TODO mhpws me endeixiInitialGeometricD...
        {
            ox_i_shell_midsurface = new double[8][];
            for (int j = 0; j < 8; j++)
            {
                ox_i_shell_midsurface[j] = new double[] { element.INodes[j].X, element.INodes[j].Y, element.INodes[j].Z, };                               
            }
        }

        

        private void UpdateCoordinateDataForDirectVectorsAndMidsurface(double[] localdisplacements, out double[][] tx_i_shell_midsurface)
        {
            double ak;
            double bk;
            tx_i_shell_midsurface = new double[8][];
            for (int j = 0; j < 8; j++)
            {
                tx_i_shell_midsurface[j] = new double[3];
            }

            for (int k = 0; k < 8; k++)
            {
                for (int l = 0; l < 3; l++)
                {
                    tx_i_shell_midsurface[k][l] = ox_i_shell_midsurface[k][l] + localdisplacements[5 * k + l];
                }
                //update twn tU kai tUvec 
                tU[k][0] = localdisplacements[5 * k + 0];
                tU[k][1] = localdisplacements[5 * k + 1];
                tU[k][2] = localdisplacements[5 * k + 2];
                ak = localdisplacements[5 * k + 3] - ak_total[k];
                ak_total[k] = localdisplacements[5 * k + 3];
                bk = localdisplacements[5 * k + 4] - bk_total[k];
                bk_total[k] = localdisplacements[5 * k + 4];
                Shell8DirectionVectorUtilities.RotateNodalDirectionVectors(ak, bk, k,  tU,  tUvec);
                // update twn tU kai tUvec ews edw                
            }
            // shmeio print dedomenwn gia debug
            ////PrintUtilities.ConvertAndWriteToFileVector(tU, @"C:\Users\turbo-x\Desktop\cohesive_check_MSOLVE_2\paradeigma_apo_arxika_swsta_shell_orthi_gia_check_tou_neou_class\orthi\CopyApoTaShellNewLoadCaseArgurhs\unused_anest_kai_t\tU_local_msolve1.txt");
            ////PrintUtilities.ConvertAndWriteToFileVector(tUvec, @"C:\Users\turbo-x\Desktop\cohesive_check_MSOLVE_2\paradeigma_apo_arxika_swsta_shell_orthi_gia_check_tou_neou_class\orthi\CopyApoTaShellNewLoadCaseArgurhs\unused_anest_kai_t\tUvec_local_msolve1.txt");
        }

        // voithitikes metavlhtes gia thn peristrofh
        //private double[] tdtVn = new double [3];
        //private double[] tdtV1 = new double[3];
        //private double[] tdtV2 = new double[3];
        //private double theta;
        //private double[] theta_vec = new double[3];
        //private double[,] s_k = new double[3, 3];        

        // prokatarktikes Methodoi kai metavlhtes apo to cohesive 16 node

        //private double[] a_12g;
        //private double a_1g;
        //private double a_2g;
        //private double[][] N1;
        //private double[][,] N3;
        //private double[][] N1_ksi;
        //private double[][] N1_heta;
        //private double[] N_i;
        //private double[] N_i_ksi;
        //private double[] N_i_heta;

        //private double ksi;
        //private double heta;
        //private double zeta;
        //private int npoint;



        //private void CalculateShapeFunctionAndGaussPointData()
        //{
        //    // PROSTHIKI RAM
        //    double a_1g=0;
        //    double a_2g=0;
        //    double[] N_i;
        //    double[] N_i_ksi;
        //    double[] N_i_heta;
        //    double ksi=0;
        //    double heta=0;
        //    //double zeta;
        //    int npoint;


        //    nGaussPoints = gp_d1_coh * gp_d2_coh;
        //    a_12g = new double[nGaussPoints];
        //    N_i = new double[8]; // 4 gia ligoterous komvous coh8
        //    N_i_ksi = new double[8]; // 
        //    N_i_heta = new double[8]; // 
        //    N1 = new double[nGaussPoints][];
        //    N3 = new double[nGaussPoints][,];
        //    N1_ksi = new double[nGaussPoints][];
        //    N1_heta = new double[nGaussPoints][];
        //    for (int l = 0; l < nGaussPoints; l++)
        //    {
        //        N1[l] = new double[8];
        //        N3[l] = new double[3, 24];
        //        N1_ksi[l] = new double[8];
        //        N1_heta[l] = new double[8];
        //    }
        //    for (int j = 0; j < gp_d1_coh; j++)
        //    {
        //        for (int k = 0; k < gp_d2_coh; k++)
        //        {
        //            npoint = j * gp_d1_coh + k;
        //            if (gp_d1_coh == 3)
        //            {
        //                ksi = 0.5 * (j - 1) * (j - 2) * (-0.774596669241483) + (-1) * (j) * (j - 2) * (0) + 0.5 * (j) * (j - 1) * (0.774596669241483);
        //                a_1g = 0.5 * (j - 1) * (j - 2) * (0.555555555555555) + (-1) * (j) * (j - 2) * (0.888888888888888) + 0.5 * (j) * (j - 1) * (0.555555555555555);
        //            }
        //            if (gp_d1_coh == 2)
        //            {
        //                ksi = (-0.577350269189626) * (j - 1) * (-1) + (0.577350269189626) * (j) * (+1);
        //                a_1g = 1;
        //            }
        //            if (gp_d2_coh == 3)
        //            {
        //                heta = 0.5 * (k - 1) * (k - 2) * (-0.774596669241483) + (-1) * (k) * (k - 2) * (0) + 0.5 * (k) * (k - 1) * (0.774596669241483);
        //                a_2g = 0.5 * (k - 1) * (k - 2) * (0.555555555555555) + (-1) * (k) * (k - 2) * (0.888888888888888) + 0.5 * (k) * (k - 1) * (0.555555555555555);
        //            }
        //            if (gp_d2_coh == 2)
        //            {
        //                heta = (-0.577350269189626) * (k - 1) * (-1) + (0.577350269189626) * (k) * (+1);
        //                a_2g = 1;
        //            }

        //            a_12g[npoint] = a_1g * a_2g;

        //            N_i[4] = 0.5 * (1 - Math.Pow(ksi, 2)) * (1 + heta);
        //            N_i[5] = 0.5 * (1 - Math.Pow(heta, 2)) * (1 - ksi);
        //            N_i[6] = 0.5 * (1 - Math.Pow(ksi, 2)) * (1 - heta);
        //            N_i[7] = 0.5 * (1 - Math.Pow(heta, 2)) * (1 + ksi);

        //            N_i[0] = 0.25 * (1 + ksi) * (1 + heta) - 0.5 * N_i[4] - 0.5 * N_i[7];
        //            N_i[1] = 0.25 * (1 - ksi) * (1 + heta) - 0.5 * N_i[4] - 0.5 * N_i[5];
        //            N_i[2] = 0.25 * (1 - ksi) * (1 - heta) - 0.5 * N_i[5] - 0.5 * N_i[6];
        //            N_i[3] = 0.25 * (1 + ksi) * (1 - heta) - 0.5 * N_i[6] - 0.5 * N_i[7];

        //            N_i_ksi[4] = (-ksi) * (1 + heta);
        //            N_i_ksi[5] = -0.5 * (1 - Math.Pow(heta, 2));
        //            N_i_ksi[6] = 0.5 * (-2 * ksi) * (1 - heta);
        //            N_i_ksi[7] = 0.5 * (1 - Math.Pow(heta, 2));
        //            N_i_ksi[0] = +0.25 * (1 + heta) - 0.5 * N_i_ksi[4] - 0.5 * N_i_ksi[7];
        //            N_i_ksi[1] = -0.25 * (1 + heta) - 0.5 * N_i_ksi[4] - 0.5 * N_i_ksi[5];
        //            N_i_ksi[2] = -0.25 * (1 - heta) - 0.5 * N_i_ksi[5] - 0.5 * N_i_ksi[6];
        //            N_i_ksi[3] = +0.25 * (1 - heta) - 0.5 * N_i_ksi[6] - 0.5 * N_i_ksi[7];

        //            N_i_heta[4] = 0.5 * (1 - Math.Pow(ksi, 2));
        //            N_i_heta[5] = 0.5 * (-2 * heta) * (1 - ksi);
        //            N_i_heta[6] = 0.5 * (1 - Math.Pow(ksi, 2)) * (-1);
        //            N_i_heta[7] = 0.5 * (-2 * heta) * (1 + ksi);
        //            N_i_heta[0] = +0.25 * (1 + ksi) - 0.5 * N_i_heta[4] - 0.5 * N_i_heta[7];
        //            N_i_heta[1] = +0.25 * (1 - ksi) - 0.5 * N_i_heta[4] - 0.5 * N_i_heta[5];
        //            N_i_heta[2] = -0.25 * (1 - ksi) - 0.5 * N_i_heta[5] - 0.5 * N_i_heta[6];
        //            N_i_heta[3] = -0.25 * (1 + ksi) - 0.5 * N_i_heta[6] - 0.5 * N_i_heta[7];

        //            for (int l = 0; l < 8; l++)  // to 8 ginetai 4 gia to cohesive8node
        //            { N1[npoint][l] = N_i[l]; }

        //            for (int l = 0; l < 3; l++)  // arxika mhdenismos twn stoixweiwn tou pinaka
        //            {
        //                for (int m = 0; m < 24; m++)
        //                { N3[npoint][l, m] = 0; }
        //            }

        //            for (int l = 0; l < 3; l++)
        //            {
        //                for (int m = 0; m < 8; m++)
        //                { N3[npoint][l, l + 3 * m] = N_i[m]; }
        //            }

        //            for (int l = 0; l < 8; l++)
        //            { N1_ksi[npoint][l] = N_i_ksi[l]; }

        //            for (int l = 0; l < 8; l++)
        //            { N1_heta[npoint][l] = N_i_heta[l]; }

        //        }
        //    }
        //    endeixiShapeFunctionAndGaussPointData = 2;
        //}

        //private double[] Get_a_12g()
        //{
        //    if (endeixiShapeFunctionAndGaussPointData == 1)
        //    {
        //        CalculateShapeFunctionAndGaussPointData();
        //        return a_12g;
        //    }
        //    else
        //    { return a_12g; }
        //}

        //private double[][] GetN1()
        //{
        //    if (endeixiShapeFunctionAndGaussPointData == 1)
        //    {
        //        CalculateShapeFunctionAndGaussPointData();
        //        return N1;
        //    }
        //    else
        //    { return N1; }
        //}

        //private double[][,] GetN3()
        //{
        //    if (endeixiShapeFunctionAndGaussPointData == 1)
        //    {
        //        CalculateShapeFunctionAndGaussPointData();
        //        return N3;
        //    }
        //    else
        //    { return N3; }
        //}

        //private double[][] GetN1_ksi()
        //{
        //    if (endeixiShapeFunctionAndGaussPointData == 1)
        //    {
        //        CalculateShapeFunctionAndGaussPointData();
        //        return N1_ksi;
        //    }
        //    else
        //    { return N1_ksi; }
        //}

        //private double[][] GetN1_heta()
        //{
        //    if (endeixiShapeFunctionAndGaussPointData == 1)
        //    {
        //        CalculateShapeFunctionAndGaussPointData();
        //        return N1_heta;
        //    }
        //    else
        //    { return N1_heta; }
        //}


        
        

        // vasikes methodoi apo to 16 node kai prosarmogh gia tous extra vathmous eleftherias

        private void GetInitialGeometricDataAndInitializeMatrices(IElement element)
        {
            //prosarmogh methodou
            (tU, tUvec) = Shell8DirectionVectorUtilities.GetInitialDirectionVectorValues(oVn_i);
            this.GetInitialGeometricDataForMidsurface(element); 
            ox_i = new double[16][];
            if ( ShellElementSide == 0)
            {
                for (int j = 0; j < 8; j++)
                { ox_i[j] = new double[] { ox_i_shell_midsurface[j][0]-0.5* tk[j] * tU[j][3], ox_i_shell_midsurface[j][1] - 0.5 * tk[j] * tU[j][4],
                                             ox_i_shell_midsurface[j][2]-0.5* tk[j] * tU[j][5], };}
                for (int j = 8; j < 16; j++)
                {ox_i[j] = new double[] { element.INodes[j].X, element.INodes[j].Y, element.INodes[j].Z, };}
            }
            else
            {
                for (int j = 0; j < 8; j++)
                {ox_i[j] = new double[] { element.INodes[j + 8].X, element.INodes[j + 8].Y, element.INodes[j + 8].Z, };}
                for (int j = 8; j < 16; j++)
                {ox_i[j] = new double[] { ox_i_shell_midsurface[j-8][0]+0.5* tk[j-8] * tU[j-8][3], ox_i_shell_midsurface[j-8][1] + 0.5 * tk[j-8] * tU[j-8][4],
                                             ox_i_shell_midsurface[j-8][2]+0.5* tk[j-8] * tU[j-8][5], };}
            }
            // prosarmogh methodou ews edw

            x_local = new double[48];
            //u_prok = new double[3, 8];
            //x_pavla = new double[3, 8];
            //k_stoixeiou_coh = new double[48, 48]; //allazei sto cohesive 8 node
            //fxk1_coh = new double[48];
            //d_trial = new double[nGaussPoints];
            //e_ksi = new double[3];
            //e_heta = new double[3];
            //e_1 = new double[3];
            //e_2 = new double[3];
            //e_3 = new double[3];
            //u = new double[3];
            nGaussPoints = gp_d1_coh * gp_d2_coh;
            //Delta = new double[nGaussPoints][];
            //R = new double[nGaussPoints][,];
            //c_1 = new double[nGaussPoints][];
            //coh_det_J_t = new double[nGaussPoints];
            //sunt_olokl = new double[nGaussPoints];
            for (int j = 0; j < nGaussPoints; j++)
            {
                //Delta[j] = new double[3];
                //R[j] = new double[3, 3];
                //c_1[j] = new double[3];
            }
            //D_tan = new double[nGaussPoints][,];
            //T_int = new double[nGaussPoints][];
            //for (int j = 0; j < nGaussPoints; j++)
            //{
            //    //D_tan[j] = new double[3, 3];
            //    //T_int[j] = new double[3];
            //}
            //M = new double[24, 24];
            //r_int = new double[48];
            //r_int_1 = new double[24];
            //RN3 = new double[nGaussPoints][,];
            //D_tan_sunt_ol = new double[3, 3];
            //D_RN3_sunt_ol = new double[3, 24];
            //T_int_sunt_ol = new double[3];
            //for (int j = 0; j < nGaussPoints; j++)
            //{
            //    RN3[j] = new double[3, 24];  //tha ginei [3,12] sto cohesive 8 node
            //}
            //this.InitializeTMatrix_and_necessary_lists_for_TKT_multiplications(element);
       }

        //private double [][] UpdateCoordinateData(double[] localdisplacements) // sto shell8disp sto calculate forces kaleitai me this.UpdateCoordinateData(localTotalDisplacements);
        //{
        //    //PROSTHIKI RAM 
        //    double[] x_local = new double[48]; // to dianusma x ths matlab sunarthshs pou einai apo t_x_global_pr
        //    double[,] u_prok = new double[3, 8];
        //    double[,] x_pavla = new double[3, 8];

        //    double[] e_ksi = new double[3];
        //    double e_ksi_norm;
        //    double[] e_heta = new double[3];
        //    double[] e_1 = new double[3];
        //    double[] e_2 = new double[3];
        //    double[] e_3 = new double[3];
        //    double e_3_norm;
        //    double[] u = new double[3];

        //    double[] coh_det_J_t = new double[nGaussPoints];

        //    double [][] Delta = new double[nGaussPoints][];
        //    double [][]c_1 = new double[nGaussPoints][];
        //    for (int j = 0; j < nGaussPoints; j++)
        //    {
        //        Delta[j] = new double[3];
        //        c_1[j] = new double[3];
        //    }

        //    double [][,] R = new double[nGaussPoints][,]; // an theloume na kanoume save to R afto tha oristhei exw sto initialize matrices opou twra to RN3
        //    for (int j = 0; j < nGaussPoints; j++)
        //    { R[j] = new double[3, 3]; }

        //    //prosarmogh methodou
        //    this.UpdateCoordinateDataForDirectVectorsAndMidsurface(localdisplacements);
        //    // ananewsi tou x_local
        //    if (endeixi_element_2==0)
        //    {
        //        for (int j = 0; j < 8; j++)
        //        {
        //            for (int k = 0; k < 3; k++)
        //            {
        //                x_local[3 * j + k] = tx_i_shell_midsurface[j][k] - 0.5 * tk[j] * tU[j][3+k];
        //            }
        //        }
        //        for (int j = 8; j < 16; j++)
        //        {
        //            for (int k = 0; k < 3; k++)
        //            {
        //                x_local[3 * j + k] = ox_i[j][k] + localdisplacements[40 + 3 * (j-8) + k];
        //            }
        //        }
        //    }
        //    else
        //    {
        //        for (int j = 0; j < 8; j++)
        //        {
        //            for (int k = 0; k < 3; k++)
        //            {
        //                x_local[3 * j + k] = ox_i[j][k] + localdisplacements[40 + 3 * j + k];
        //            }
        //        }
        //        for (int j = 8; j < 16; j++)
        //        {
        //            for (int k = 0; k < 3; k++)
        //            {
        //                x_local[3 * j + k] = tx_i_shell_midsurface[j - 8][k] + 0.5 * tk[j - 8] * tU[j - 8][3 + k];
        //            }
        //        }
        //    }
        //    // prosarmogh methodou ews edw
        //    for (int j = 0; j < 8; j++)
        //    {
        //        for (int k = 0; k < 3; k++)
        //        {
        //            u_prok[k, j] = x_local[k + 3 * j] - x_local[24 + k + 3 * j];
        //            x_pavla[k, j] = x_local[k + 3 * j] + x_local[24 + k + 3 * j];
        //        }
        //    }
        //    // sunexeia ews upologismou tou Delta gia ola ta gp

        //    for (int npoint1 = 0; npoint1 < nGaussPoints; npoint1++)
        //    {

        //        for (int l = 0; l < 3; l++)
        //        {
        //            e_ksi[l] = 0;
        //            e_heta[l] = 0;
        //            for (int m = 0; m < 8; m++) // tha ginei 4 sto cohesive 8 node
        //            {
        //                e_ksi[l] += N1_ksi[npoint1][m] * x_pavla[l, m];
        //                e_heta[l] += N1_heta[npoint1][m] * x_pavla[l, m];
        //            }
        //            e_ksi[l] = 0.5 * e_ksi[l];
        //            e_heta[l] = 0.5 * e_heta[l];
        //        }
        //        this.cross(e_ksi, e_heta, e_3);
        //        e_3_norm = Math.Sqrt(e_3[0] * e_3[0] + e_3[1] * e_3[1] + e_3[2] * e_3[2]);
        //        e_ksi_norm = Math.Sqrt(e_ksi[0] * e_ksi[0] + e_ksi[1] * e_ksi[1] + e_ksi[2] * e_ksi[2]);
        //        for (int l = 0; l < 3; l++)
        //        {
        //            e_3[l] = e_3[l] / e_3_norm;
        //            e_1[l] = e_ksi[l] / e_ksi_norm;
        //        }
        //        this.cross(e_1, e_3, e_2);
        //        for (int l = 0; l < 3; l++)
        //        {
        //            R[npoint1][l, 0] = e_1[l];
        //            R[npoint1][l, 1] = e_2[l];
        //            R[npoint1][l, 2] = e_3[l];

        //        }
        //        for (int l = 0; l < 3; l++)
        //        { u[l] = 0; }
        //        for (int l = 0; l < 3; l++)
        //        {
        //            for (int m = 0; m < 8; m++)  // pithanws gia to cohesive 8
        //            {
        //                u[l] += u_prok[l, m] * N1[npoint1][m];
        //            }
        //        }
        //        //for (int l = 0; l < 3; l++)
        //        //{ Delta[npoint1][l] = 0; }    //sp1
        //        for (int l = 0; l < 3; l++)
        //        {
        //            for (int m = 0; m < 3; m++)
        //            {
        //                Delta[npoint1][l] += R[npoint1][m, l] * u[m];
        //            }
        //        }

        //        this.cross(e_ksi, e_heta, c_1[npoint1]);
        //        coh_det_J_t[npoint1] = Math.Sqrt(c_1[npoint1][0] * c_1[npoint1][0] + c_1[npoint1][1] * c_1[npoint1][1] + c_1[npoint1][2] * c_1[npoint1][2]);
        //        sunt_olokl[npoint1] = coh_det_J_t[npoint1] * a_12g[npoint1];

                
        //        // upologimsos tou RN3 edw anti entos tou initializeRN3 kai tou updateForces
        //        RN3[npoint1] = new double [3,24];                
        //        for (int l = 0; l < 3; l++)
        //        {
        //            for (int m = 0; m < 24; m++)
        //            {
        //                for (int n = 0; n < 3; n++)
        //                { RN3[npoint1][l, m] += R[npoint1][l, n] * N3[npoint1][n, m]; }
        //            }
        //        }



        //    }
        //    //this.UpdateTMatrix();
        //    print_counter += 1;           
        //    return Delta;
        //}

        private double[][] UpdateCoordinateDataAndCalculateDisplacementVector(double[] localdisplacements) // sto shell8disp sto calculate forces kaleitai me this.UpdateCoordinateData(localTotalDisplacements);
        {
            // retrieve twn metavlhtwn pou den tha apothikevontai (shape function data)
            //(double[][] N1old, double[][,] N3old, double[][] N1_ksi, double[][] N1_heta, double[] a_12g) =interpolation.GetShapeFunctionAndGaussPointData(nGaussPoints, gp_d1_coh, gp_d2_coh);
            IReadOnlyList<Matrix2D> shapeFunctionDerivatives = interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForStiffness);
            IReadOnlyList<Vector> N1 = interpolation.EvaluateFunctionsAtGaussPoints(QuadratureForStiffness);
            IReadOnlyList<Matrix2D> N3 = interpolation.EvaluateN3ShapeFunctionsReorganized(QuadratureForStiffness);

            //PROSTHIKI RAM 
            //double[] x_local = new double[48]; // to dianusma x ths matlab sunarthshs pou einai apo t_x_global_pr
            double[,] u_prok = new double[3, 8];
            double[,] x_pavla = new double[3, 8];

            double[] e_ksi = new double[3];
            double e_ksi_norm;
            double[] e_heta = new double[3];
            double[] e_1 = new double[3];
            double[] e_2 = new double[3];
            double[] e_3 = new double[3];
            double e_3_norm;
            double[] u = new double[3];

            double[] coh_det_J_t = new double[nGaussPoints];

            double[][] Delta = new double[nGaussPoints][];
            double[][] c_1 = new double[nGaussPoints][];
            for (int j = 0; j < nGaussPoints; j++)
            {
                Delta[j] = new double[3];
                c_1[j] = new double[3];
            }

            double[][,] R = new double[nGaussPoints][,]; // an theloume na kanoume save to R afto tha oristhei exw sto initialize matrices opou twra to RN3
            for (int j = 0; j < nGaussPoints; j++)
            { R[j] = new double[3, 3]; }

            //prosarmogh methodou
            this.UpdateCoordinateDataForDirectVectorsAndMidsurface(localdisplacements,out double [][] tx_i_shell_midsurface);
            // ananewsi tou x_local
            if (ShellElementSide == 0)
            {
                for (int j = 0; j < 8; j++)
                {
                    for (int k = 0; k < 3; k++)
                    {
                        x_local[3 * j + k] = tx_i_shell_midsurface[j][k] - 0.5 * tk[j] * tU[j][3 + k];
                    }
                }
                for (int j = 8; j < 16; j++)
                {
                    for (int k = 0; k < 3; k++)
                    {
                        x_local[3 * j + k] = ox_i[j][k] + localdisplacements[40 + 3 * (j - 8) + k];
                    }
                }
            }
            else
            {
                for (int j = 0; j < 8; j++)
                {
                    for (int k = 0; k < 3; k++)
                    {
                        x_local[3 * j + k] = ox_i[j][k] + localdisplacements[40 + 3 * j + k];
                    }
                }
                for (int j = 8; j < 16; j++)
                {
                    for (int k = 0; k < 3; k++)
                    {
                        x_local[3 * j + k] = tx_i_shell_midsurface[j - 8][k] + 0.5 * tk[j - 8] * tU[j - 8][3 + k];
                    }
                }
            }
            // prosarmogh methodou ews edw
            for (int j = 0; j < 8; j++)
            {
                for (int k = 0; k < 3; k++)
                {
                    u_prok[k, j] = x_local[k + 3 * j] - x_local[24 + k + 3 * j];
                    x_pavla[k, j] = x_local[k + 3 * j] + x_local[24 + k + 3 * j];
                }
            }
            // sunexeia ews upologismou tou Delta gia ola ta gp

            for (int npoint1 = 0; npoint1 < nGaussPoints; npoint1++)
            {

                for (int l = 0; l < 3; l++)
                {
                    e_ksi[l] = 0;
                    e_heta[l] = 0;
                    for (int m = 0; m < 8; m++) // tha ginei 4 sto cohesive 8 node
                    {
                        e_ksi[l] += shapeFunctionDerivatives[npoint1][0, m] * x_pavla[l, m];
                        e_heta[l] += shapeFunctionDerivatives[npoint1][1, m] * x_pavla[l, m];
                    }
                    e_ksi[l] = 0.5 * e_ksi[l];
                    e_heta[l] = 0.5 * e_heta[l];
                }
                this.Cross(e_ksi, e_heta, e_3);
                e_3_norm = Math.Sqrt(e_3[0] * e_3[0] + e_3[1] * e_3[1] + e_3[2] * e_3[2]);
                e_ksi_norm = Math.Sqrt(e_ksi[0] * e_ksi[0] + e_ksi[1] * e_ksi[1] + e_ksi[2] * e_ksi[2]);
                for (int l = 0; l < 3; l++)
                {
                    e_3[l] = e_3[l] / e_3_norm;
                    e_1[l] = e_ksi[l] / e_ksi_norm;
                }
                this.Cross(e_1, e_3, e_2);
                for (int l = 0; l < 3; l++)
                {
                    R[npoint1][l, 0] = e_1[l];
                    R[npoint1][l, 1] = e_2[l];
                    R[npoint1][l, 2] = e_3[l];

                }
                for (int l = 0; l < 3; l++)
                { u[l] = 0; }
                for (int l = 0; l < 3; l++)
                {
                    for (int m = 0; m < 8; m++)  // pithanws gia to cohesive 8
                    {
                        u[l] += u_prok[l, m] * N1[npoint1][m];
                    }
                }
                //for (int l = 0; l < 3; l++)
                //{ Delta[npoint1][l] = 0; }    //sp1
                for (int l = 0; l < 3; l++)
                {
                    for (int m = 0; m < 3; m++)
                    {
                        Delta[npoint1][l] += R[npoint1][m, l] * u[m];
                    }
                }

                //this.cross(e_ksi, e_heta, c_1[npoint1]);
                //coh_det_J_t[npoint1] = Math.Sqrt(c_1[npoint1][0] * c_1[npoint1][0] + c_1[npoint1][1] * c_1[npoint1][1] + c_1[npoint1][2] * c_1[npoint1][2]);
                //sunt_olokl[npoint1] = coh_det_J_t[npoint1] * a_12g[npoint1];


                //// upologimsos tou RN3 edw anti entos tou initializeRN3 kai tou updateForces
                //RN3[npoint1] = new double[3, 24];
                //for (int l = 0; l < 3; l++)
                //{
                //    for (int m = 0; m < 24; m++)
                //    {
                //        for (int n = 0; n < 3; n++)
                //        { RN3[npoint1][l, m] += R[npoint1][l, n] * N3[npoint1][n, m]; }
                //    }
                //}



            }
            //this.UpdateTMatrix();
            //print_counter += 1;
            return Delta;
        }

        private Tuple<Matrix2D[], double[]> CalculateNecessaryMatricesForStiffnessMatrixAndForcesVectorCalculations() // sto shell8disp sto calculate forces kaleitai me this.UpdateCoordinateData(localTotalDisplacements);
        {
            // retrieve twn metavlhtwn pou den tha apothikevontai (shape function data)
            //(double[][] N1old, double[][,] N3old, double[][] N1_ksi, double[][] N1_heta, double[] a_12g) = interpolation.GetShapeFunctionAndGaussPointData(nGaussPoints, gp_d1_coh, gp_d2_coh);

            IReadOnlyList<Vector> N1 = interpolation.EvaluateFunctionsAtGaussPoints(QuadratureForStiffness);
            IReadOnlyList<Matrix2D> N3 = interpolation.EvaluateN3ShapeFunctionsReorganized(QuadratureForStiffness);
            IReadOnlyList<Matrix2D> shapeFunctionDerivatives = interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForStiffness);

            //initialize twn mhtrwwn pou ta epistrefontai
            double[]sunt_olokl = new double[nGaussPoints];
            Matrix2D []  RtN3 = new Matrix2D [nGaussPoints];

            //PROSTHIKI RAM 
            //double[] x_local = new double[48]; // to dianusma x ths matlab sunarthshs pou einai apo t_x_global_pr
            //double[,] u_prok = new double[3, 8];
            double[,] x_pavla = new double[3, 8];

            
            double[] e_1 = new double[3];
            double[] e_2 = new double[3];
            double[] e_3 = new double[3];
            double e_3_norm;
            //double[] u = new double[3];

            double[] coh_det_J_t = new double[nGaussPoints];

            //double[][] Delta = new double[nGaussPoints][];
            double[][] c_1 = new double[nGaussPoints][];
            for (int j = 0; j < nGaussPoints; j++)
            {
                //Delta[j] = new double[3];
                c_1[j] = new double[3];
            }

            Matrix2D[] R = new Matrix2D [nGaussPoints]; // an theloume na kanoume save to R afto tha oristhei exw sto initialize matrices opou twra to RN3
            for (int j = 0; j < nGaussPoints; j++)
            { R[j] = new Matrix2D( 3, 3); }

            //prosarmogh methodou
            //this.UpdateCoordinateDataForDirectVectorsAndMidsurface(localdisplacements);
            // ananewsi tou x_local
            //if (endeixi_element_2 == 0)
            //{
            //    for (int j = 0; j < 8; j++)
            //    {
            //        for (int k = 0; k < 3; k++)
            //        {
            //            x_local[3 * j + k] = tx_i_shell_midsurface[j][k] - 0.5 * tk[j] * tU[j][3 + k];
            //        }
            //    }
            //    for (int j = 8; j < 16; j++)
            //    {
            //        for (int k = 0; k < 3; k++)
            //        {
            //            x_local[3 * j + k] = ox_i[j][k] + localdisplacements[40 + 3 * (j - 8) + k];
            //        }
            //    }
            //}
            //else
            //{
            //    for (int j = 0; j < 8; j++)
            //    {
            //        for (int k = 0; k < 3; k++)
            //        {
            //            x_local[3 * j + k] = ox_i[j][k] + localdisplacements[40 + 3 * j + k];
            //        }
            //    }
            //    for (int j = 8; j < 16; j++)
            //    {
            //        for (int k = 0; k < 3; k++)
            //        {
            //            x_local[3 * j + k] = tx_i_shell_midsurface[j - 8][k] + 0.5 * tk[j - 8] * tU[j - 8][3 + k];
            //        }
            //    }
            //}
            // prosarmogh methodou ews edw
            for (int j = 0; j < 8; j++)
            {
                for (int k = 0; k < 3; k++)
                {
                    //u_prok[k, j] = x_local[k + 3 * j] - x_local[24 + k + 3 * j]; //den xreiazetai sth metodo nesMAtrcForStiffnAndForce
                    x_pavla[k, j] = x_local[k + 3 * j] + x_local[24 + k + 3 * j];
                }
            }
            // sunexeia ews upologismou tou Delta gia ola ta gp

            for (int npoint1 = 0; npoint1 < nGaussPoints; npoint1++)
            {
                double[] e_ksi = new double[3];                
                double[] e_heta = new double[3];
                for (int l = 0; l < 3; l++)
                {
                    //e_ksi[l] = 0;
                    //e_heta[l] = 0;
                    for (int m = 0; m < 8; m++) // tha ginei 4 sto cohesive 8 node
                    {
                        e_ksi[l] += shapeFunctionDerivatives[npoint1][0, m] * x_pavla[l, m];
                        e_heta[l] += shapeFunctionDerivatives[npoint1][1, m] * x_pavla[l, m];
                    }
                    e_ksi[l] = 0.5 * e_ksi[l];
                    e_heta[l] = 0.5 * e_heta[l];
                }
                this.Cross(e_ksi, e_heta, e_3);
                e_3_norm = Math.Sqrt(e_3[0] * e_3[0] + e_3[1] * e_3[1] + e_3[2] * e_3[2]);
                double e_ksi_norm = Math.Sqrt(e_ksi[0] * e_ksi[0] + e_ksi[1] * e_ksi[1] + e_ksi[2] * e_ksi[2]);
                for (int l = 0; l < 3; l++)
                {
                    e_3[l] = e_3[l] / e_3_norm;
                    e_1[l] = e_ksi[l] / e_ksi_norm;
                }
                this.Cross(e_1, e_3, e_2);
                for (int l = 0; l < 3; l++)
                {
                    R[npoint1][l, 0] = e_1[l];
                    R[npoint1][l, 1] = e_2[l];
                    R[npoint1][l, 2] = e_3[l];

                }
                //for (int l = 0; l < 3; l++)
                //{ u[l] = 0; }
                //for (int l = 0; l < 3; l++)
                //{
                //    for (int m = 0; m < 8; m++)  // pithanws gia to cohesive 8
                //    {
                //        u[l] += u_prok[l, m] * N1[npoint1][m];
                //    }
                //}
                ////for (int l = 0; l < 3; l++)
                ////{ Delta[npoint1][l] = 0; }    //sp1
                //for (int l = 0; l < 3; l++)
                //{
                //    for (int m = 0; m < 3; m++)
                //    {
                //        Delta[npoint1][l] += R[npoint1][m, l] * u[m];
                //    }
                //}

                this.Cross(e_ksi, e_heta, c_1[npoint1]);
                coh_det_J_t[npoint1] = Math.Sqrt(c_1[npoint1][0] * c_1[npoint1][0] + c_1[npoint1][1] * c_1[npoint1][1] + c_1[npoint1][2] * c_1[npoint1][2]);
                sunt_olokl[npoint1] = coh_det_J_t[npoint1] * QuadratureForStiffness.IntegrationPoints[npoint1].Weight;


                // upologimsos tou RN3 edw anti entos tou initializeRN3 kai tou updateForces

                RtN3[npoint1] = R[npoint1].Transpose() * N3[npoint1];
                //RtN3[npoint1] = new double[3, 24];
                //for (int l = 0; l < 3; l++)
                //{
                //    for (int m = 0; m < 24; m++)
                //    {
                //        for (int n = 0; n < 3; n++)
                //        { RtN3[npoint1][l, m] += R[npoint1][n,l] * N3[npoint1][n, m]; }
                //    }
                //}



            }
            //this.UpdateTMatrix();
            //print_counter += 1;
            Tuple<Matrix2D[], double[]> RtN3AndSunt_olokl = new Tuple<Matrix2D[], double[]>(RtN3, sunt_olokl);
            return RtN3AndSunt_olokl;
        }




        private void Cross(double[] A, double[] B, double[] C)
        {
            C[0] = A[1] * B[2] - A[2] * B[1];
            C[1] = A[2] * B[0] - A[0] * B[2];
            C[2] = A[0] * B[1] - A[1] * B[0];
        }


        // methodos dhmiourgia tou T to opoio tha xrhsimopoihthei gia thn TKT kai to opoio tha einai etoimo gia xrhsh gia 
        //to prwto efaptomeniko mhtrwo duskampsias

        //methodos ananewsh toy T h opoia tha kaleitai amesws meta apo thn update coordinate data h mporoume na thn kaloume
        // sthn teleftaia grammh mesa sthn update coordinate data 

        // methodoi CalcualteForces{ lipsi dedomenwn apo material-Updateforces-return fxk1_coh}
        // kai Stifness Matrix {lipsi dedomenwn apo material- UpdateKmatrixes- return IMatrix element_stifnessMatrix}
        // prepei se kapoia apo tis parapanw diadikasies na mpei h praxh TKT kai TF

        // diadikasia T_int= .... pou ginetai mesa sto calcForces kai sto stifnessMatrix h antistoixh den htan allh telika apo
        //th lhpsh apo to material twn aparaithtwn dedomenwn

        //pollaplasiasmos mhtrwnn programatismos edw h mporei kai oxi


        //private double[,] T;
        //private double[,] eye3;

        //private double[,] Kii_A;
        //private double[,] k_stoixeiou_coh2;
        //private double[] fxk2_coh;

        //private void InitializeTMatrix_and_necessary_lists_for_TKT_multiplications(Element element)
        //{
        //    //Kii_A = new double[24, 40];
        //    //k_stoixeiou_coh2 = new double[64, 64];
        //    //fxk2_coh = new double[64];

        //    T = new double[24, 40];
        //    double[,] eye3 = new double[3, 3];
        //    for (int m = 0; m < 3; m++)
        //    {
        //        for (int l = 0; l < 3; l++)
        //        {
        //            eye3[m, l] = 0;
        //        }
        //        eye3[m, m] = 1;

        //    }
        //    for (int m = 0; m < 8; m++)
        //    {
        //        for (int n = 0; n < 3; n++)
        //        {
        //            for (int l = 0; l < 3; l++)
        //            {
        //                T[3 * m + n, 5 * m + l] = eye3[n, l];
        //            }
        //        }
                
        //    }
        //    if(endeixi_element_2==0)
        //    {
        //        for (int m = 0; m < 8; m++)
        //        {
        //            for (int n = 0; n < 3; n++)
        //            {
        //                T[3 * m + n, 5 * m + 3] =0.5*tk[m]*tUvec[m][3 + n];
        //                T[3 * m + n, 5 * m + 4] = -0.5 * tk[m] * tUvec[m][n];
        //            }
        //        }
        //    }
        //    else
        //    {
        //        for (int m = 0; m < 8; m++)
        //        {
        //            for (int n = 0; n < 3; n++)
        //            {
        //                T[3 * m + n, 5 * m + 3] = -0.5 * tk[m] * tUvec[m][3 + n];
        //                T[3 * m + n, 5 * m + 4] = +0.5 * tk[m] * tUvec[m][n];
        //            }
        //        }
        //    }

        //}

        private double [,] CalculateTMatrix(IElement element)
        {
            //Kii_A = new double[24, 40];
            //k_stoixeiou_coh2 = new double[64, 64];
            //fxk2_coh = new double[64];

            double[,] T;
            T = new double[24, 40];
            double[,] eye3 = new double[3, 3];
            for (int m = 0; m < 3; m++)
            {                
                eye3[m, m] = 1;
            }
            for (int m = 0; m < 8; m++)
            {
                for (int n = 0; n < 3; n++)
                {
                    for (int l = 0; l < 3; l++)
                    {
                        T[3 * m + n, 5 * m + l] = eye3[n, l];
                    }
                }

            }
            if (ShellElementSide == 0)
            {
                for (int m = 0; m < 8; m++)
                {
                    for (int n = 0; n < 3; n++)
                    {
                        T[3 * m + n, 5 * m + 3] = 0.5 * tk[m] * tUvec[m][3 + n];
                        T[3 * m + n, 5 * m + 4] = -0.5 * tk[m] * tUvec[m][n];
                    }
                }
            }
            else
            {
                for (int m = 0; m < 8; m++)
                {
                    for (int n = 0; n < 3; n++)
                    {
                        T[3 * m + n, 5 * m + 3] = -0.5 * tk[m] * tUvec[m][3 + n];
                        T[3 * m + n, 5 * m + 4] = +0.5 * tk[m] * tUvec[m][n];
                    }
                }
            }

            return T;
        }

        //private void UpdateTMatrix()
        //{
        //    if (endeixi_element_2 == 0)
        //    {
        //        for (int m = 0; m < 8; m++)
        //        {
        //            for (int n = 0; n < 3; n++)
        //            {
        //                T[3 * m + n, 5 * m + 3] = 0.5 * tk[m] * tUvec[m][3 + n];
        //                T[3 * m + n, 5 * m + 4] = -0.5 * tk[m] * tUvec[m][n];
        //            }
        //        }
        //    }
        //    else
        //    {
        //        for (int m = 0; m < 8; m++)
        //        {
        //            for (int n = 0; n < 3; n++)
        //            {
        //                T[3 * m + n, 5 * m + 3] = -0.5 * tk[m] * tUvec[m][3 + n];
        //                T[3 * m + n, 5 * m + 4] = +0.5 * tk[m] * tUvec[m][n];
        //            }
        //        }
        //    }
        //    //
        //    //if (print_counter == 2)
        //    //{
        //    //    PrintUtilities.WriteToFile(T, @"C:\Users\turbo-x\Desktop\cohesive_check_MSOLVE_2\paradeigma_apo_arxika_swsta_shell_orthi_gia_check_tou_neou_class\orthi\CopyApoTaShellNewLoadCaseArgurhs\unused_anest_kai_t\T_output_1.txt");
        //    //}
        //    //if (print_counter == 3)
        //    //{
        //    //    PrintUtilities.WriteToFile(T, @"C:\Users\turbo-x\Desktop\cohesive_check_MSOLVE_2\paradeigma_apo_arxika_swsta_shell_orthi_gia_check_tou_neou_class\orthi\CopyApoTaShellNewLoadCaseArgurhs\unused_anest_kai_t\T_output_2.txt");
        //    //}
        //    //
        //}

        // telikh morfh pinakwn gia embeding kai endiameses metavlhtes


        private double [] MultiplyForcesForEmbedding(double [] fxk1_coh,Element element)
        {
            //upologismos tou T
            double[,] T;
            T = CalculateTMatrix(element);
            // declare tou dianusmatos pou tha epistrefei h methodos
            double [] fxk2_coh = new double[64];
            if (ShellElementSide==0)
            {
                for (int n = 0; n < 40; n++)
                {
                    //fxk2_coh[n] = 0; //sp1
                    for (int p = 0; p < 24; p++)
                    {
                        fxk2_coh[n] += T[p, n] * fxk1_coh[p];
                    }
                }
                for (int n = 40; n <64 ; n++)
                {
                    fxk2_coh[n] = fxk1_coh[n - 16]; // ta oria tou loop einai swsta exoume 24 epanalipseis omws kai to n-16 swsto (px to prwto 40-16=24 sto fxk1)
                }

            }
            else
            {
                for (int n = 0; n < 24; n++)
                {
                    //fxk2_coh[n] = fxk1_coh[n];
                    fxk2_coh[40+n] = fxk1_coh[n];
                }
                for (int n = 0; n < 40; n++)
                {
                    //fxk2_coh[n+24] = 0;
                    //fxk2_coh[n] = 0;  //sp1 afto einai sp1 , to apo panw einai lathos
                    for (int p = 0; p < 24; p++)
                    {
                        //fxk2_coh[n+24] += T[p, n] * fxk1_coh[p+24];
                        fxk2_coh[n] += T[p, n] * fxk1_coh[p + 24];
                    }
                }
            }

            //PrintUtilities.WriteToFile(T,
            //       @"C:\Users\turbo-x\Desktop\cohesive_check_MSOLVE_2\paradeigma_apo_arxika_swsta_shell_orthi_gia_check_tou_neou_class\orthi\CopyApoTaShellNewLoadCaseArgurhs\unused_anest_kai_t\T_updated_A_1__mh_an.txt");
            return fxk2_coh;
        }

        private double [,] MultiplyStifnessMatrixForEmbedding(double [,] k_stoixeiou_coh, IElement element)
        {
            //upologismos tou T
            double[,] T;
            T = CalculateTMatrix( element);
            // declare tou mhtrwou pou tha epistrefei h methodos kathws kai enos aparaithtou gia tous upologismous
            double [,] k_stoixeiou_coh2 = new double[64, 64];
            double [,] Kii_A = new double[24, 40];
            if (ShellElementSide == 0)
            {
                //upologismos Kii_A (mhdenismos kai upologismoi)
                //for (int n = 0; n < 24; n++)
                //{
                //    for (int p = 0; p < 40; p++)
                //    {
                //        Kii_A[n, p] = 0;
                //    }
                //} //sp1
                for (int n = 0; n < 24; n++)
                {
                    for (int p = 0; p < 40; p++)
                    {
                        for (int k = 0; k < 24; k++)
                        {
                            Kii_A[n, p] += k_stoixeiou_coh[n, k] * T[k, p];
                        }
                    }
                }
                // upologismos perioxhs 11 tou Tt_K_T (mhdenismos kai upologismoi)
                //for (int n = 0; n < 40; n++)
                //{
                //    for (int p = 0; p < 40; p++)
                //    {
                //        k_stoixeiou_coh2[n, p] = 0;
                //    }
                //} //sp1
                for (int n = 0; n < 40; n++)
                {
                    for (int p = 0; p < 40; p++)
                    {
                        for (int k = 0; k < 24; k++)
                        {
                            k_stoixeiou_coh2[n, p] += T[k, n] * Kii_A[k, p];
                        }
                    }
                }
                // upologismos perioxhs 12 tou Tt_K_T (mhdenismos kai upologismoi)
                //for (int n = 0; n < 40; n++)
                //{
                //    for (int p = 0; p < 24; p++)
                //    {
                //        k_stoixeiou_coh2[n, 40 + p] = 0;
                //    }
                //}//sp1
                for (int n = 0; n < 40; n++)
                {
                    for (int p = 0; p < 24; p++)
                    {
                        for (int k = 0; k < 24; k++)
                        {
                            k_stoixeiou_coh2[n, 40 + p] += T[k, n] * k_stoixeiou_coh[k, 24 + p];
                        }
                    }
                }
                // upologismos perioxhs 21 tou Tt_K_T (mhdenismos kai upologismoi)
                //for (int n = 0; n < 24; n++)
                //{
                //    for (int p = 0; p < 40; p++)
                //    {
                //        k_stoixeiou_coh2[40+n, p] = 0;
                //    }
                //} //sp1
                for (int n = 0; n < 24; n++)
                {
                    for (int p = 0; p < 40; p++)
                    {
                        for (int k = 0; k < 24; k++)
                        {
                            k_stoixeiou_coh2[40+n, p] +=  k_stoixeiou_coh[24+n, k] * T[k, p];
                        }
                    }
                }
                // upologismos perioxhs 22 tou Tt_K_T (copy paste apo k_stoixeiou_coh)
                for (int n = 0; n < 24; n++)
                {
                    for (int p = 0; p < 24; p++)
                    {
                        k_stoixeiou_coh2[40 + n,40 + p] = k_stoixeiou_coh[24+n,24+p];
                    }
                }
                //PrintUtilities.SeparateAndWriteToFile(k_stoixeiou_coh2,
                //    @"C:\Users\turbo-x\Desktop\cohesive_check_MSOLVE_2\paradeigma_apo_arxika_swsta_shell_orthi_gia_check_tou_neou_class\orthi\CopyApoTaShellNewLoadCaseArgurhs\unused_anest_kai_t\K_stoixeiou_coh_2_A_1.txt",
                //    @"C:\Users\turbo-x\Desktop\cohesive_check_MSOLVE_2\paradeigma_apo_arxika_swsta_shell_orthi_gia_check_tou_neou_class\orthi\CopyApoTaShellNewLoadCaseArgurhs\unused_anest_kai_t\K_stoixeiou_coh_2_B_1.txt");
                //PrintUtilities.WriteToFile(k_stoixeiou_coh_2_A, @"C:\Users\turbo-x\Desktop\cohesive_check_MSOLVE_2\paradeigma_apo_arxika_swsta_shell_orthi_gia_check_tou_neou_class\orthi\CopyApoTaShellNewLoadCaseArgurhs\unused_anest_kai_t\K_stoixeiou_coh_2_A_1.txt");
                //PrintUtilities.WriteToFile(k_stoixeiou_coh_2_B, @"C:\Users\turbo-x\Desktop\cohesive_check_MSOLVE_2\paradeigma_apo_arxika_swsta_shell_orthi_gia_check_tou_neou_class\orthi\CopyApoTaShellNewLoadCaseArgurhs\unused_anest_kai_t\K_stoixeiou_coh_2_B_1.txt");
                //PrintUtilities.WriteToFile(k_stoixeiou_coh2, @"C:\Users\turbo-x\Desktop\cohesive_check_MSOLVE_2\apotelesmata_MSOLVE\abc_1.txt");
                //PrintUtilities.WriteToFile(k_stoixeiou_coh2, @"C:\Users\turbo-x\Desktop\cohesive_check_MSOLVE_2\paradeigma_apo_arxika_swsta_shell_orthi_gia_check_tou_neou_class\orthi\CopyApoTaShellNewLoadCaseArgurhs\unused_anest_kai_t\abc_1.txt");
            }
            else
            {
                //upologismos Kii_A pou en prokeimenw einai K22_A (mhdenismos kai upologismoi)
                //for (int n = 0; n < 24; n++)
                //{
                //    for (int p = 0; p < 40; p++)
                //    {
                //        Kii_A[n, p] = 0;
                //    }
                //}
                for (int n = 0; n < 24; n++)
                {
                    for (int p = 0; p < 40; p++)
                    {
                        for (int k = 0; k < 24; k++)
                        {
                            Kii_A[n, p] += k_stoixeiou_coh[24+n, 24+k] * T[k, p];
                        }
                    }
                }
                // upologismos perioxhs 11 tou Tt_K_T (copy paste apo k_stoixeiou_coh) -->22
                for (int n = 0; n < 24; n++)
                {
                    for (int p = 0; p < 24; p++)
                    {
                        //k_stoixeiou_coh2[n,p] = k_stoixeiou_coh[n,p];
                        k_stoixeiou_coh2[40+n,40+ p] = k_stoixeiou_coh[n, p];
                    }
                }
                // upologismos perioxhs 12 tou Tt_K_T (mhdenismos kai upologismoi) --> 21
                //for (int n = 0; n < 24; n++)
                //{
                //    for (int p = 0; p < 40; p++)
                //    {
                //        //k_stoixeiou_coh2[ n, 24+p] = 0;
                //        k_stoixeiou_coh2[40+n, p] = 0;
                //    }
                //} //sp1
                for (int n = 0; n < 24; n++)
                {
                    for (int p = 0; p < 40; p++)
                    {
                        for (int k = 0; k < 24; k++)
                        {
                            //k_stoixeiou_coh2[n,24+ p] += k_stoixeiou_coh[n, 24+k] * T[k, p];
                            k_stoixeiou_coh2[40+n, p] += k_stoixeiou_coh[n, 24 + k] * T[k, p];
                        }
                    }
                }
                // upologismos perioxhs 21 tou Tt_K_T (mhdenismos kai upologismoi) -->12
                //for (int n = 0; n < 40; n++)
                //{
                //    for (int p = 0; p < 24; p++)
                //    {
                //        //k_stoixeiou_coh2[24+n, p] = 0;
                //        k_stoixeiou_coh2[n,40+ p] = 0;
                //    }
                //}
                for (int n = 0; n < 40; n++)
                {
                    for (int p = 0; p < 24; p++)
                    {
                        for (int k = 0; k < 24; k++)
                        {
                            //k_stoixeiou_coh2[24+n, p] += T[k, n] * k_stoixeiou_coh[24+k, p];
                            k_stoixeiou_coh2[n,40+ p] += T[k, n] * k_stoixeiou_coh[24 + k, p];
                        }
                    }
                }
                // upologismos perioxhs 22 tou Tt_K_T (mhdenismos kai upologismoi) -->11
                //for (int n = 0; n < 40; n++)
                //{
                //    for (int p = 0; p < 40; p++)
                //    {
                //        //k_stoixeiou_coh2[24+n, 24+p] = 0;
                //        k_stoixeiou_coh2[n, p] = 0;
                //    }
                //}
                for (int n = 0; n < 40; n++)
                {
                    for (int p = 0; p < 40; p++)
                    {
                        for (int k = 0; k < 24; k++)
                        {
                            //k_stoixeiou_coh2[24 + n, 24 + p] += T[k, n] * Kii_A[k, p];
                            k_stoixeiou_coh2[ n,  p] += T[k, n] * Kii_A[k, p];
                        }
                    }
                }
                
            }
            //if (print_counter == 0)
            //{
            //    PrintUtilities.SeparateAndWriteToFile(k_stoixeiou_coh,
            //       @"C:\Users\turbo-x\Desktop\cohesive_check_MSOLVE_2\paradeigma_apo_arxika_swsta_shell_orthi_gia_check_tou_neou_class\orthi\CopyApoTaShellNewLoadCaseArgurhs\unused_anest_kai_t\K_stoixeiou_coh_arxiko_A_1_mh_an.txt",
            //       @"C:\Users\turbo-x\Desktop\cohesive_check_MSOLVE_2\paradeigma_apo_arxika_swsta_shell_orthi_gia_check_tou_neou_class\orthi\CopyApoTaShellNewLoadCaseArgurhs\unused_anest_kai_t\K_stoixeiou_coh_arxiko_B_1_mh_an.txt");
            //}
            ////if (print_counter == 1)
            ////{
            ////    PrintUtilities.SeparateAndWriteToFile(k_stoixeiou_coh,
            ////        @"C:\Users\turbo-x\Desktop\cohesive_check_MSOLVE_2\paradeigma_apo_arxika_swsta_shell_orthi_gia_check_tou_neou_class\orthi\CopyApoTaShellNewLoadCaseArgurhs\unused_anest_kai_t\K_stoixeiou_coh_updated_A_1_mh_an.txt",
            ////        @"C:\Users\turbo-x\Desktop\cohesive_check_MSOLVE_2\paradeigma_apo_arxika_swsta_shell_orthi_gia_check_tou_neou_class\orthi\CopyApoTaShellNewLoadCaseArgurhs\unused_anest_kai_t\K_stoixeiou_coh_updated_B_1_mh_an.txt");
            ////}
            //if (print_counter == 2)
            //{
            //    PrintUtilities.SeparateAndWriteToFile(k_stoixeiou_coh,
            //        @"C:\Users\turbo-x\Desktop\cohesive_check_MSOLVE_2\paradeigma_apo_arxika_swsta_shell_orthi_gia_check_tou_neou_class\orthi\CopyApoTaShellNewLoadCaseArgurhs\unused_anest_kai_t\K_stoixeiou_coh_local_updated_swsto_A_mh_an.txt",
            //        @"C:\Users\turbo-x\Desktop\cohesive_check_MSOLVE_2\paradeigma_apo_arxika_swsta_shell_orthi_gia_check_tou_neou_class\orthi\CopyApoTaShellNewLoadCaseArgurhs\unused_anest_kai_t\K_stoixeiou_coh_local_updated_swsto_B_mh_an.txt");
            //    PrintUtilities.WriteToFile(T,
            //       @"C:\Users\turbo-x\Desktop\cohesive_check_MSOLVE_2\paradeigma_apo_arxika_swsta_shell_orthi_gia_check_tou_neou_class\orthi\CopyApoTaShellNewLoadCaseArgurhs\unused_anest_kai_t\T_updated_swsto_mh_an.txt");
            //    PrintUtilities.WriteToFileVector(fxk1_coh, @"C:\Users\turbo-x\Desktop\cohesive_check_MSOLVE_2\paradeigma_apo_arxika_swsta_shell_orthi_gia_check_tou_neou_class\orthi\CopyApoTaShellNewLoadCaseArgurhs\unused_anest_kai_t\fxk1_updated_swsto_mh_an.txt");
            //}

            return k_stoixeiou_coh2;
        }

        private double[] UpdateForces(Element element, Matrix2D[] RtN3, double[] sunt_olokl)
        {
            double [] fxk2_coh = new double[64];

            //for (int j = 0; j < 48; j++) // allagh sto cohesive 8 node
            //{
            //    fxk1_coh[j] = 0;
            //}
            double [] fxk1_coh = new double[48];

            for (int npoint1 = 0; npoint1 < nGaussPoints; npoint1++)
            {
                // o upologismos tou RN3 ginetai entos tou updatecoordinateData
                //for (int l = 0; l < 3; l++)
                //{
                //    for (int m = 0; m < 24; m++)
                //    {
                //        RN3[npoint1][l, m] = 0;
                //    }
                //}
                //for (int l = 0; l < 3; l++)
                //{
                //    for (int m = 0; m < 24; m++)
                //    {
                //        for (int n = 0; n < 3; n++)
                //        { RN3[npoint1][l, m] += R[npoint1][l, n] * N3[npoint1][n, m]; }
                //    }
                //}
                double [] T_int_sunt_ol = new double[3];
                for (int l = 0; l < 3; l++)
                {
                    T_int_sunt_ol[l] = materialsAtGaussPoints[npoint1].Tractions[l] * sunt_olokl[npoint1]; //T_int[npoint1][l] * sunt_olokl[npoint1];
                }

                //for (int l = 0; l < 24; l++)
                //{ r_int_1[l] = 0; }
                double [] r_int_1 = new double[24];
                for (int l = 0; l < 24; l++)
                {
                    for (int m = 0; m < 3; m++)
                    { r_int_1[l] += RtN3[npoint1][m, l] * T_int_sunt_ol[m]; }
                }
                for (int l = 0; l < 24; l++)
                {
                    fxk1_coh[l] += r_int_1[l];
                    fxk1_coh[24 + l] += (-r_int_1[l]);
                }
            }
            //
            //if (print_counter == 2)
            //{
            //    PrintUtilities.WriteToFileVector(fxk1_coh, @"C:\Users\turbo-x\Desktop\cohesive_check_MSOLVE_2\paradeigma_apo_arxika_swsta_shell_orthi_gia_check_tou_neou_class\orthi\CopyApoTaShellNewLoadCaseArgurhs\unused_anest_kai_t\fxk1_coh_output_1.txt");
            //}
            //if (print_counter == 3)
            //{
            //    PrintUtilities.WriteToFileVector(fxk1_coh, @"C:\Users\turbo-x\Desktop\cohesive_check_MSOLVE_2\paradeigma_apo_arxika_swsta_shell_orthi_gia_check_tou_neou_class\orthi\CopyApoTaShellNewLoadCaseArgurhs\unused_anest_kai_t\fxk1_coh_output_2.txt");
            //}
            //
            fxk2_coh=this.MultiplyForcesForEmbedding(fxk1_coh,element);
            //PrintUtilities.WriteToFileVector(fxk1_coh, @"C:\Users\turbo-x\Desktop\cohesive_check_MSOLVE_2\paradeigma_apo_arxika_swsta_shell_orthi_gia_check_tou_neou_class\orthi\CopyApoTaShellNewLoadCaseArgurhs\unused_anest_kai_t\fxk1_mh_an.txt");
            ////
            //if (print_counter == 1)
            //{
            //    PrintUtilities.WriteToFileVector(fxk2_coh, @"C:\Users\turbo-x\Desktop\cohesive_check_MSOLVE_2\paradeigma_apo_arxika_swsta_shell_orthi_gia_check_tou_neou_class\orthi\CopyApoTaShellNewLoadCaseArgurhs\unused_anest_kai_t\fxk2_coh_output_0_mh_an.txt");
            //}
            //if (print_counter == 2)
            //{
            //    PrintUtilities.WriteToFileVector(fxk2_coh, @"C:\Users\turbo-x\Desktop\cohesive_check_MSOLVE_2\paradeigma_apo_arxika_swsta_shell_orthi_gia_check_tou_neou_class\orthi\CopyApoTaShellNewLoadCaseArgurhs\unused_anest_kai_t\fxk2_coh_output_1_mh_an.txt");
            //}
            //if (print_counter == 3)
            //{
            //    PrintUtilities.WriteToFileVector(fxk2_coh, @"C:\Users\turbo-x\Desktop\cohesive_check_MSOLVE_2\paradeigma_apo_arxika_swsta_shell_orthi_gia_check_tou_neou_class\orthi\CopyApoTaShellNewLoadCaseArgurhs\unused_anest_kai_t\fxk2_coh_output_2_mh_an.txt");
            //}
            //
            return fxk2_coh;
        }

        private double [,] UpdateKmatrices(IElement element, Matrix2D[] RtN3, double[] sunt_olokl)
        {
            double [,] k_stoixeiou_coh2 = new double[64, 64];
            //for (int k = 0; k < 48; k++) // allagh sto cohesive 8 node
            //{
            //    for (int j = 0; j < 48; j++)
            //    {
            //        k_stoixeiou_coh[k, j] = 0;
            //    }
            //}
            double [,] k_stoixeiou_coh = new double[48, 48];


            for (int npoint1 = 0; npoint1 < nGaussPoints; npoint1++)
            {
                Matrix2D D_tan_sunt_ol = new Matrix2D (3, 3);
                for (int l = 0; l < 3; l++)
                {
                    for (int m = 0; m < 3; m++)
                    {
                        D_tan_sunt_ol[l, m] = materialsAtGaussPoints[npoint1].ConstitutiveMatrix[l, m] * sunt_olokl[npoint1];// D_tan[npoint1][l, m] * sunt_olokl[npoint1];
                    }
                }

                //for (int l = 0; l < 3; l++)
                //{
                //    for (int m = 0; m < 24; m++)  // omoiws to 24 gia ta cohesive 8 node
                //    {
                //        D_RN3_sunt_ol[l, m] = 0;
                //    }
                //}

                Matrix2D  D_RtN3_sunt_ol = D_tan_sunt_ol*RtN3[npoint1];
                //for (int l = 0; l < 3; l++)
                //{
                //    for (int m = 0; m < 24; m++)
                //    {
                //        for (int n = 0; n < 3; n++)
                //        { D_RtN3_sunt_ol[l, m] += D_tan_sunt_ol[l, n] * RtN3[npoint1][n, m]; }
                //    }
                //}

                //for (int l = 0; l < 24; l++)
                //{
                //    for (int m = 0; m < 24; m++)
                //    {
                //        M[l, m] = 0;
                //    }
                //}

                Matrix2D M = RtN3[npoint1].Transpose() * D_RtN3_sunt_ol;
                //double [,] M = new double[24, 24];
                //for (int l = 0; l < 24; l++)
                //{
                //    for (int m = 0; m < 24; m++)
                //    {
                //        for (int n = 0; n < 3; n++)
                //        {
                //            M[l, m] += RtN3[npoint1][n, l] * D_RtN3_sunt_ol[n, m];
                //        }
                //    }
                //}

                for (int l = 0; l < 24; l++)
                {
                    for (int m = 0; m < 24; m++)
                    {
                        k_stoixeiou_coh[l, m] += M[l, m];
                        k_stoixeiou_coh[l, 24 + m] += -M[l, m];
                        k_stoixeiou_coh[24 + l, m] += -M[l, m];
                        k_stoixeiou_coh[24 + l, 24 + m] += M[l, m];
                    }
                }
            }

            k_stoixeiou_coh2 =this.MultiplyStifnessMatrixForEmbedding(k_stoixeiou_coh,element);
            //if (print_counter == 1)
            //{
            //    PrintUtilities.SeparateAndWriteToFile(k_stoixeiou_coh2,
            //            @"C:\Users\turbo-x\Desktop\cohesive_check_MSOLVE_2\paradeigma_apo_arxika_swsta_shell_orthi_gia_check_tou_neou_class\orthi\CopyApoTaShellNewLoadCaseArgurhs\unused_anest_kai_t\K_stoixeiou_coh_2_A_1_mh_an.txt",
            //            @"C:\Users\turbo-x\Desktop\cohesive_check_MSOLVE_2\paradeigma_apo_arxika_swsta_shell_orthi_gia_check_tou_neou_class\orthi\CopyApoTaShellNewLoadCaseArgurhs\unused_anest_kai_t\K_stoixeiou_coh_2_B_1_mh_an.txt");
            //}
            //if (print_counter == 2)
            //{
            //    PrintUtilities.SeparateAndWriteToFile(k_stoixeiou_coh2,
            //            @"C:\Users\turbo-x\Desktop\cohesive_check_MSOLVE_2\paradeigma_apo_arxika_swsta_shell_orthi_gia_check_tou_neou_class\orthi\CopyApoTaShellNewLoadCaseArgurhs\unused_anest_kai_t\K_stoixeiou_coh2_updated_swsto_A_mh_an.txt",
            //            @"C:\Users\turbo-x\Desktop\cohesive_check_MSOLVE_2\paradeigma_apo_arxika_swsta_shell_orthi_gia_check_tou_neou_class\orthi\CopyApoTaShellNewLoadCaseArgurhs\unused_anest_kai_t\K_stoixeiou_coh2_updated_swsto_B_mh_an.txt");
            //}
            return k_stoixeiou_coh2;
        }

        //int print_counter = 0;
        public Tuple<double[], double[]> CalculateStresses(Element element, double[] localTotalDisplacementsSuperElement, double[] localdDisplacementsSuperElement)
        {
            double[][] Delta = new double[nGaussPoints][];
            double[] localTotalDisplacements = dofEnumerator.GetTransformedDisplacementsVector(localTotalDisplacementsSuperElement); // PROSTHIKI EMBEDDED: metonomasthkan ta localTotalDisplacements se "SuperElement" apo ta opoia tha upologizontai pleon ta localTotalDisplacemeents
            //double[] localDisplacements = dofEnumerator.GetTransformedVector(localdDisplacementsSuperElement); // PROSTHIKI EMBEDDED: omoiws apla einai commented out epeidh den ta Xrhsimopoioume pouthena
            //if (print_counter == 3)
            //{
            //    //PrintUtilities.WriteToFileVector(localTotalDisplacements, @"C:\Users\turbo-x\Desktop\cohesive_check_MSOLVE_2\paradeigma_apo_arxika_swsta_shell_orthi_gia_check_tou_neou_class\orthi\CopyApoTaShellNewLoadCaseArgurhs\unused_anest_kai_t\abc_vec.txt");
            //}
            Delta=this.UpdateCoordinateDataAndCalculateDisplacementVector(localTotalDisplacements);
            for (int i = 0; i < materialsAtGaussPoints.Length; i++)
            {
                materialsAtGaussPoints[i].UpdateMaterial(Delta[i]);
            }
            return new Tuple<double[], double[]>(Delta[materialsAtGaussPoints.Length - 1], materialsAtGaussPoints[materialsAtGaussPoints.Length - 1].Tractions);
            //TODO mono to teleftaio dianusma tha epistrefei?
            
        }

        public double[] CalculateForces(Element element, double[] localTotalDisplacementsSuperElement, double[] localdDisplacementsSuperelement)
        {
            double[] fxk2_coh = new double[64];
            //double[] localTotalDisplacements = dofEnumerator.GetTransformedVector(localTotalDisplacementsSuperElement); // PROSTHIKI EMBEDDED: metonomasthkan ta localTotalDisplacements se "SuperElement" apo ta opoia tha upologizontai pleon ta localTotalDisplacemeents
            //double[] localDisplacements = dofEnumerator.GetTransformedVector(localdDisplacementsSuperElement); // PROSTHIKI EMBEDDED: omoiws apla einai commented out epeidh den ta Xrhsimopoioume pouthena
            //for (int i = 0; i < materialsAtGaussPoints.Length; i++)
            //{
            //    for (int j = 0; j < 3; j++)
            //    { T_int[i][j] = materialsAtGaussPoints[i].Stresses[j]; }
            //}
            Tuple<Matrix2D[], double[]> RtN3AndSunt_olokl;
            RtN3AndSunt_olokl = CalculateNecessaryMatricesForStiffnessMatrixAndForcesVectorCalculations();
            Matrix2D[] RtN3;
            RtN3 = RtN3AndSunt_olokl.Item1;
            double[] sunt_olokl; //osa kai ta gpoints
            sunt_olokl = RtN3AndSunt_olokl.Item2;

            fxk2_coh = this.UpdateForces(element,RtN3,sunt_olokl); //<-- sto Update afto exei enswmatwthei to Tt*fxk1_coh
            //return fxk1_coh;

            //return fxk2_coh; // PROSTHIKI EMBEDDED: Line commented out 
            return dofEnumerator.GetTransformedForcesVector(fxk2_coh);// PROSTHIKI EMBEDDED
        }

        public double[] CalculateForcesForLogging(Element element, double[] localDisplacements)
        {
            return CalculateForces(element, localDisplacements, new double[localDisplacements.Length]);
        }

        public virtual IMatrix2D StiffnessMatrix(IElement element)   //anti IMatrix2D<double> kai sta upoloipa compa
        {
            double [,] k_stoixeiou_coh2 = new double[64, 64];
            if (MatrixIsNotInitialized)
            {
                //this.CalculateShapeFunctionAndGaussPointData();
                this.GetInitialGeometricDataAndInitializeMatrices(element);
                this.UpdateCoordinateDataAndCalculateDisplacementVector(new double[64]); // mporei na lamvanetai apo to update... to return tou (pou einai Delta[][] p.x. gia material)
                MatrixIsNotInitialized = false;
            }
            //for (int i = 0; i < materialsAtGaussPoints.Length; i++)
            //{
            //    for (int j = 0; j < 3; j++)
            //    {
            //        for (int k = 0; k < 3; k++)
            //        { D_tan[i][j, k] = materialsAtGaussPoints[i].ConstitutiveMatrix[j, k]; }
            //    }
            //}
            Tuple<Matrix2D[], double[]> RtN3AndSunt_olokl;
            RtN3AndSunt_olokl = CalculateNecessaryMatricesForStiffnessMatrixAndForcesVectorCalculations();
            Matrix2D[] RtN3;
            RtN3 = RtN3AndSunt_olokl.Item1;
            double[] sunt_olokl; //osa kai ta gpoints
            sunt_olokl = RtN3AndSunt_olokl.Item2;

            k_stoixeiou_coh2 =this.UpdateKmatrices(element, RtN3, sunt_olokl); //<-- sto Update afto exei mpei to Tt*K*T 
            //IMatrix2D<double> element_stiffnessMatrix = new Matrix2D<double>(k_stoixeiou_coh); // TODO giati de ginetai return dof.Enumerator.GetTransformedMatrix, xrhsh symmetric
            IMatrix2D element_stiffnessMatrix = new Matrix2D(k_stoixeiou_coh2); // TODO giati de ginetai return dof.Enumerator.GetTransformedMatrix, xrhsh symmetric
            //return element_stiffnessMatrix; //PROSTHIKI EMBEDDED: Line Commented Out
            IMatrix2D superElement_stiffnessMatrix= dofEnumerator.GetTransformedMatrix(element_stiffnessMatrix);// PROSTHIKI GIA CHECK EMBEDDDED
            return dofEnumerator.GetTransformedMatrix(element_stiffnessMatrix); // PROSTHIKI EMBEDDED
        }


        public bool MaterialModified
        {
            get
            {
                foreach (ICohesiveZoneMaterial3D material in materialsAtGaussPoints)
                    if (material.Modified) return true;
                return false;
            }
        }

        public void ResetMaterialModified()
        {
            foreach (ICohesiveZoneMaterial3D material in materialsAtGaussPoints) material.ResetModified();
        }

        public void ClearMaterialState()
        {
            foreach (ICohesiveZoneMaterial3D m in materialsAtGaussPoints) m.ClearState();
        }

        public void SaveMaterialState()
        {
            foreach (ICohesiveZoneMaterial3D m in materialsAtGaussPoints) m.SaveState();
        }

        public void ClearMaterialStresses()
        {
            foreach (ICohesiveZoneMaterial3D m in materialsAtGaussPoints) m.ClearTractions();
        }

        // omoiws me hexa 8 shell8disp implemented
        public int ID
        {
            get { return 13; }
        }
        public ElementDimensions ElementDimensions
        {
            get { return ElementDimensions.ThreeD; }
        }

        public IElementDOFEnumerator DOFEnumerator
        {
            get { return dofEnumerator; }
            set { dofEnumerator = value; }
        }

        public virtual IList<IList<DOFType>> GetElementDOFTypes(IElement element)
        {
            return dofTypes;
        }

        // aplopoihtika implemented mhdenikes masses gia cohesive - not implemented
        public double[] CalculateAccelerationForces(Element element, IList<MassAccelerationLoad> loads)
        {
            return new double[64];
        }

        public virtual IMatrix2D MassMatrix(IElement element)
        {
            return new Matrix2D(64, 64);
        }

        public virtual IMatrix2D DampingMatrix(IElement element)
        {

            return new Matrix2D(64, 64);
        }

        #region EMBEDDED
        private readonly List<EmbeddedNode> embeddedNodes = new List<EmbeddedNode>(); //
        public IList<EmbeddedNode> EmbeddedNodes { get { return embeddedNodes; } } // opws Beam3D

        public Dictionary<DOFType, int> GetInternalNodalDOFs(Element element, Node node)//
        {
            int index = 0;
            foreach (var elementNode in element.Nodes)
            {
                if (node.ID == elementNode.ID)
                    break;
                index++;
            }
            if (index >= 16)
                throw new ArgumentException(String.Format("GetInternalNodalDOFs: Node {0} not found in element {1}.", node.ID, element.ID));

            if (index >= 8)
            {
                int index2 = index - 8;
                return new Dictionary<DOFType, int>() { { DOFType.X, 39+3*index2+1 }, { DOFType.Y, 39 + 3 * index2 + 2 }, { DOFType.Z, 39 + 3 * index2 + 3 } };
            }
            else
            {
                return new Dictionary<DOFType, int>() { { DOFType.X, + 5 * index + 0 }, { DOFType.Y, + 5 * index + 1 }, { DOFType.Z, + 5 * index + 2 },
                                                        { DOFType.RotX, + 5 * index + 3 }, { DOFType.RotY, + 5 * index + 4 }};
            }          
        }

        public double[] GetLocalDOFValues(Element hostElement, double[] hostDOFValues) // omoiws Beam3D
        {
            //if (transformation == null)
            //    throw new InvalidOperationException("Requested embedded node values for element that has no embedded nodes.");
            //if (hostElementList == null)
            //    throw new InvalidOperationException("Requested host element list for element that has no embedded nodes.");
            //int index = hostElementList.IndexOf(hostElement);
            //if (index < 0)
            //    throw new ArgumentException("Requested host element is not inside host element list.");

            //double[] values = new double[transformation.Columns];
            //int multiplier = hostElement.ElementType.DOFEnumerator.GetDOFTypes(hostElement).SelectMany(d => d).Count();
            //int vectorIndex = 0;
            //for (int i = 0; i < index; i++)
            //    vectorIndex += isNodeEmbedded[i] ? 3 : multiplier;
            //Array.Copy(hostDOFValues, 0, values, vectorIndex, multiplier);

            //return (transformation * new Vector<double>(values)).Data;

            return dofEnumerator.GetTransformedDisplacementsVector(hostDOFValues);
        }
        #endregion
    }
}
