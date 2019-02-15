using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.MultiscaleAnalysis.SupportiveClasses
{
    public class renumbering
    {
        public int[] sunol_nodes_numbering { get; set; }

        public renumbering()
        {

        }
        public renumbering(int[] sunol_nodes_numbering)
        {
            this.sunol_nodes_numbering = sunol_nodes_numbering;
        }

        public int GetNewNodeNumbering(int initial_node_number)
        {
            return sunol_nodes_numbering[initial_node_number - 1];
        }
    }

    public class rveMatrixParameters
    {
        public double E_disp { get; set; }
        public double ni_disp { get; set; }
        public double L01 { get; set; }
        public double L02 { get; set; }
        public double L03 { get; set; }
        public int hexa1 { get; set; }
        public int hexa2 { get; set; }
        public int hexa3 { get; set; }

        public rveMatrixParameters()
        {

        }
        public rveMatrixParameters(double E_disp, double ni_disp, double L01, double L02, double L03, int hexa1, int hexa2, int hexa3)
        {
            this.E_disp = E_disp;
            this.ni_disp = ni_disp;
            this.L01 = L01;
            this.L02 = L02;
            this.L03 = L03;
            this.hexa1 = hexa1;
            this.hexa2 = hexa2;
            this.hexa3 = hexa3;
        }
    }

    public class o_x_parameters
    {
        //public double E_disp { get; set; }
        //public double ni_disp { get; set; }
        //public double L01 { get; set; }
        //public double L02 { get; set; }
        //public double L03 { get; set; }
        //public int hexa1 { get; set; }
        //public int hexa2 { get; set; }
        //public int hexa3 { get; set; }

        public o_x_parameters()
        {

        }
        public o_x_parameters(double E_disp, double ni_disp, double L01, double L02, double L03, int hexa1, int hexa2, int hexa3)
        {
            //this.E_disp = E_disp;
            //this.ni_disp = ni_disp;
            //this.L01 = L01;
            //this.L02 = L02;
            //this.L03 = L03;
            //this.hexa1 = hexa1;
            //this.hexa2 = hexa2;
            //this.hexa3 = hexa3;
        }
    }

    public class grapheneSheetParameters
    {
        // parametroi shell
        public double E_shell; // GPa = 1000Mpa = 1000N / mm2
        public double ni_shell; // stathera poisson
        public int elem1;
        public int elem2;
        public double L1;// nm
        public double L2;// nm
        public double L3; // nm
        public double a1_shell; // nm
        public double tk;  // 0.0125016478913782nm
                           //parametroi cohesive epifaneias
        public double T_o_3;// Gpa = 1000Mpa = 1000N / mm2
        public double D_o_3; // nm
        public double D_f_3; // nm
        public double T_o_1;// Gpa
        public double D_o_1; // nm
        public double D_f_1; // nm
        public double n_curve = 1.4;

        public grapheneSheetParameters()
        {

        }
        public grapheneSheetParameters(double E_shell, double ni_shell, int elem1, int elem2, double L1, double L2, double L3, double a1_shell, double tk,
            double T_o_3, double D_o_3, double D_f_3, double T_o_1, double D_o_1, double D_f_1, double n_curve)
        {
            this.E_shell = E_shell; // GPa = 1000Mpa = 1000N / mm2
            this.ni_shell = ni_shell; // stathera poisson
            this.elem1 = elem1;
            this.elem2 = elem2;
            this.L1 = L1;// nm
            this.L2 = L2;// nm
            this.L3 = L3; // nm
            this.a1_shell = a1_shell; // nm
            this.tk = tk;  // 0.0125016478913782nm

            //parametroi cohesive epifaneias
            //T_o_3, D_o_3,D_f_3,T_o_1,D_o_1,D_f_1,n_curve
            this.T_o_3 = T_o_3;// Gpa = 1000Mpa = 1000N / mm2
            this.D_o_3 = D_o_3; // nm
            this.D_f_3 = D_f_3; // nm

            this.T_o_1 = T_o_1;// Gpa
            this.D_o_1 = D_o_1; // nm
            this.D_f_1 = D_f_1; // nm

            this.n_curve = n_curve;
        }
    }
}
