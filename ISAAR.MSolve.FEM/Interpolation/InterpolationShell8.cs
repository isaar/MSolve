using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.FEM.Interpolation
{
    public class InterpolationShell8
    {
        private static readonly InterpolationShell8 uniqueInstance = new InterpolationShell8();

        public static InterpolationShell8 UniqueInstance => uniqueInstance;
        
        public (double[][] gaussCoordinates, double[][] shapeFunctions, double[][] shapeFunctionDerivatives, double[] a_123g) 
            GetShapeFunctions(int gp_d1, int gp_d2, int gp_d3)
        {
            double[][] gaussCoordinates;//3 dianysmata me tis timew tvn ksi heta zeta se ola ta gauss points
            double[][] shapeFunctions;// 8 dianusmata me tis times twn N1....N8 se kathe gauss point
            double[][] shapeFunctionDerivatives;// 16 dianusmata me tis times twn N1ksi....N8ksi,N1heta,....N8heta se kathe gauss point

            double ksi = 0;
            double heta = 0;
            double zeta = 0;
            double a_1g = 0;
            double a_2g = 0;
            double a_3g = 0;

            int nGaussPoints = gp_d1 * gp_d2 * gp_d3;
            double[] a_123g = new double[nGaussPoints];
            gaussCoordinates = new double[3][];
            for (int l = 0; l < 3; l++)
            { gaussCoordinates[l] = new double[nGaussPoints]; }
            for (int l = 0; l < gp_d3; l++)
            {
                for (int k = 0; k < gp_d2; k++)
                {
                    for (int j = 0; j < gp_d1; j++)
                    {
                        int npoint = l * (gp_d1 * gp_d2) + k * gp_d1 + j;
                        if (gp_d1 == 3)
                        {
                            ksi = 0.5 * (j - 1) * (j - 2) * (-0.774596669241483) + (-1) * (j) * (j - 2) * (0) + 0.5 * (j) * (j - 1) * (0.774596669241483);
                            a_1g = 0.5 * (j - 1) * (j - 2) * (0.555555555555555) + (-1) * (j) * (j - 2) * (0.888888888888888) + 0.5 * (j) * (j - 1) * (0.555555555555555);
                        }
                        if (gp_d1 == 2)
                        {
                            ksi = (-0.577350269189626) * (j - 1) * (-1) + (0.577350269189626) * (j) * (+1);
                            a_1g = 1;
                        }
                        if (gp_d2 == 3)
                        {
                            heta = 0.5 * (k - 1) * (k - 2) * (-0.774596669241483) + (-1) * (k) * (k - 2) * (0) + 0.5 * (k) * (k - 1) * (0.774596669241483);
                            a_2g = 0.5 * (k - 1) * (k - 2) * (0.555555555555555) + (-1) * (k) * (k - 2) * (0.888888888888888) + 0.5 * (k) * (k - 1) * (0.555555555555555);
                        }
                        if (gp_d2 == 2)
                        {
                            heta = (-0.577350269189626) * (k - 1) * (-1) + (0.577350269189626) * (k) * (+1);
                            a_2g = 1;
                        }
                        if (gp_d3 == 3)
                        {
                            zeta = 0.5 * (l - 1) * (l - 2) * (-0.774596669241483) + (-1) * (l) * (l - 2) * (0) + 0.5 * (l) * (l - 1) * (0.774596669241483);
                            a_3g = 0.5 * (l - 1) * (l - 2) * (0.555555555555555) + (-1) * (l) * (l - 2) * (0.888888888888888) + 0.5 * (l) * (l - 1) * (0.555555555555555);
                        }
                        if (gp_d3 == 2)
                        {
                            zeta = (-0.577350269189626) * (l - 1) * (-1) + (0.577350269189626) * (l) * (+1);
                            a_3g = 1;
                        }
                        gaussCoordinates[0][npoint] = ksi;
                        gaussCoordinates[1][npoint] = heta;
                        gaussCoordinates[2][npoint] = zeta;

                        a_123g[npoint] = a_1g * a_2g * a_3g;
                    }
                }
            }

            shapeFunctions = new double[8][];
            for (int j = 0; j < 8; j++)
            { shapeFunctions[j] = new double[nGaussPoints]; }
            for (int j = 0; j < nGaussPoints; j++)
            {
                shapeFunctions[4][j] = 0.5 * (1 - Math.Pow(gaussCoordinates[0][j], 2)) * (1 + gaussCoordinates[1][j]);
                shapeFunctions[5][j] = 0.5 * (1 - Math.Pow(gaussCoordinates[1][j], 2)) * (1 - gaussCoordinates[0][j]);
                shapeFunctions[6][j] = 0.5 * (1 - Math.Pow(gaussCoordinates[0][j], 2)) * (1 - gaussCoordinates[1][j]);
                shapeFunctions[7][j] = 0.5 * (1 - Math.Pow(gaussCoordinates[1][j], 2)) * (1 + gaussCoordinates[0][j]);
                shapeFunctions[0][j] = 0.25 * (1 + gaussCoordinates[0][j]) * (1 + gaussCoordinates[1][j]) - 0.5 * shapeFunctions[4][j] - 0.5 * shapeFunctions[7][j];
                shapeFunctions[1][j] = 0.25 * (1 - gaussCoordinates[0][j]) * (1 + gaussCoordinates[1][j]) - 0.5 * shapeFunctions[4][j] - 0.5 * shapeFunctions[5][j];
                shapeFunctions[2][j] = 0.25 * (1 - gaussCoordinates[0][j]) * (1 - gaussCoordinates[1][j]) - 0.5 * shapeFunctions[5][j] - 0.5 * shapeFunctions[6][j];
                shapeFunctions[3][j] = 0.25 * (1 + gaussCoordinates[0][j]) * (1 - gaussCoordinates[1][j]) - 0.5 * shapeFunctions[6][j] - 0.5 * shapeFunctions[7][j];
            }



            shapeFunctionDerivatives = new double[16][];
            for (int j = 0; j < 16; j++)
            { shapeFunctionDerivatives[j] = new double[nGaussPoints]; }
            for (int j = 0; j < nGaussPoints; j++)
            {
                //Ni_ksi
                shapeFunctionDerivatives[4][j] = (-gaussCoordinates[0][j]) * (1 + gaussCoordinates[1][j]);
                shapeFunctionDerivatives[5][j] = -0.5 * (1 - Math.Pow(gaussCoordinates[1][j], 2));
                shapeFunctionDerivatives[6][j] = 0.5 * (-2 * gaussCoordinates[0][j]) * (1 - gaussCoordinates[1][j]);
                shapeFunctionDerivatives[7][j] = 0.5 * (1 - Math.Pow(gaussCoordinates[1][j], 2));
                shapeFunctionDerivatives[0][j] = +0.25 * (1 + gaussCoordinates[1][j]) - 0.5 * shapeFunctionDerivatives[4][j] - 0.5 * shapeFunctionDerivatives[7][j];
                shapeFunctionDerivatives[1][j] = -0.25 * (1 + gaussCoordinates[1][j]) - 0.5 * shapeFunctionDerivatives[4][j] - 0.5 * shapeFunctionDerivatives[5][j];
                shapeFunctionDerivatives[2][j] = -0.25 * (1 - gaussCoordinates[1][j]) - 0.5 * shapeFunctionDerivatives[5][j] - 0.5 * shapeFunctionDerivatives[6][j];
                shapeFunctionDerivatives[3][j] = +0.25 * (1 - gaussCoordinates[1][j]) - 0.5 * shapeFunctionDerivatives[6][j] - 0.5 * shapeFunctionDerivatives[7][j];
                //Ni_heta
                shapeFunctionDerivatives[12][j] = 0.5 * (1 - Math.Pow(gaussCoordinates[0][j], 2));
                shapeFunctionDerivatives[13][j] = 0.5 * (-2 * gaussCoordinates[1][j]) * (1 - gaussCoordinates[0][j]);
                shapeFunctionDerivatives[14][j] = 0.5 * (1 - Math.Pow(gaussCoordinates[0][j], 2)) * (-1);
                shapeFunctionDerivatives[15][j] = 0.5 * (-2 * gaussCoordinates[1][j]) * (1 + gaussCoordinates[0][j]);
                shapeFunctionDerivatives[8][j] = +0.25 * (1 + gaussCoordinates[0][j]) - 0.5 * shapeFunctionDerivatives[12][j] - 0.5 * shapeFunctionDerivatives[15][j];
                shapeFunctionDerivatives[9][j] = +0.25 * (1 - gaussCoordinates[0][j]) - 0.5 * shapeFunctionDerivatives[12][j] - 0.5 * shapeFunctionDerivatives[13][j];
                shapeFunctionDerivatives[10][j] = -0.25 * (1 - gaussCoordinates[0][j]) - 0.5 * shapeFunctionDerivatives[13][j] - 0.5 * shapeFunctionDerivatives[14][j];
                shapeFunctionDerivatives[11][j] = -0.25 * (1 + gaussCoordinates[0][j]) - 0.5 * shapeFunctionDerivatives[14][j] - 0.5 * shapeFunctionDerivatives[15][j];
            }
           
            return (gaussCoordinates, shapeFunctions, shapeFunctionDerivatives, a_123g);
        }

    }
}
