using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.FEM.Interpolation
{
    public class InterpolationHexa8Reverse
    {
        private static readonly InterpolationHexa8Reverse uniqueInstance = new InterpolationHexa8Reverse();

        public static InterpolationHexa8Reverse UniqueInstance => uniqueInstance;

        public (double[,] Ni_ksi, double[,] Ni_heta, double[,] Ni_zeta, double[] a_123g) GetShapeFunctionDerivatives(int gp_d1_disp, int gp_d2_disp, int gp_d3_disp)
        {
            double[] a_123g;
            double a_1g = 0;
            double a_2g = 0;
            double a_3g = 0;
            double ksi = 0;
            double heta = 0;
            double zeta = 0;
            int npoint;
            //double[,] Ni;
            double[,] Ni_ksi;
            double[,] Ni_heta;
            double[,] Ni_zeta;

            int nGaussPoints = gp_d1_disp * gp_d2_disp * gp_d3_disp;
            a_123g = new double[nGaussPoints];

            //Ni = new double[8, nGaussPoints]; // den sxetizetai me ta coh elements alla
            Ni_ksi = new double[8, nGaussPoints]; // me to prokat_disp (einai shapefunctionData)
            Ni_heta = new double[8, nGaussPoints]; // 
            Ni_zeta = new double[8, nGaussPoints];

            for (int l = 0; l < gp_d3_disp; l++)
            {
                for (int k = 0; k < gp_d2_disp; k++)
                {
                    for (int j = 0; j < gp_d1_disp; j++)
                    {
                        npoint = l * (gp_d1_disp * gp_d2_disp) + k * gp_d1_disp + j;
                        if (gp_d1_disp == 3)
                        {
                            ksi = 0.5 * (j - 1) * (j - 2) * (-0.774596669241483) + (-1) * (j) * (j - 2) * (0) + 0.5 * (j) * (j - 1) * (0.774596669241483);
                            a_1g = 0.5 * (j - 1) * (j - 2) * (0.555555555555555) + (-1) * (j) * (j - 2) * (0.888888888888888) + 0.5 * (j) * (j - 1) * (0.555555555555555);
                        }
                        if (gp_d1_disp == 2)
                        {
                            ksi = (-0.577350269189626) * (j - 1) * (-1) + (0.577350269189626) * (j) * (+1);
                            a_1g = 1;
                        }
                        if (gp_d2_disp == 3)
                        {
                            heta = 0.5 * (k - 1) * (k - 2) * (-0.774596669241483) + (-1) * (k) * (k - 2) * (0) + 0.5 * (k) * (k - 1) * (0.774596669241483);
                            a_2g = 0.5 * (k - 1) * (k - 2) * (0.555555555555555) + (-1) * (k) * (k - 2) * (0.888888888888888) + 0.5 * (k) * (k - 1) * (0.555555555555555);
                        }
                        if (gp_d2_disp == 2)
                        {
                            heta = (-0.577350269189626) * (k - 1) * (-1) + (0.577350269189626) * (k) * (+1);
                            a_2g = 1;
                        }
                        if (gp_d3_disp == 3)
                        {
                            zeta = 0.5 * (l - 1) * (l - 2) * (-0.774596669241483) + (-1) * (l) * (l - 2) * (0) + 0.5 * (l) * (l - 1) * (0.774596669241483);
                            a_3g = 0.5 * (l - 1) * (l - 2) * (0.555555555555555) + (-1) * (l) * (l - 2) * (0.888888888888888) + 0.5 * (l) * (l - 1) * (0.555555555555555);
                        }
                        if (gp_d3_disp == 2)
                        {
                            zeta = (-0.577350269189626) * (l - 1) * (-1) + (0.577350269189626) * (l) * (+1);
                            a_3g = 1;
                        }
                        a_123g[npoint] = a_1g * a_2g * a_3g;

                        //Ni[0, npoint] = 0.125 * (1 + ksi) * (1 + heta) * (1 + zeta); //N1(ksi,heta,zeta) sto sugkekrimeno GP "npoint" sto [,]
                        //Ni[1, npoint] = 0.125 * (1 - ksi) * (1 + heta) * (1 + zeta); //N2(ksi,heta,zeta) sto sugkekrimeno GP "npoint" sto [,]
                        //Ni[2, npoint] = 0.125 * (1 - ksi) * (1 - heta) * (1 + zeta);
                        //Ni[3, npoint] = 0.125 * (1 + ksi) * (1 - heta) * (1 + zeta);
                        //Ni[4, npoint] = 0.125 * (1 + ksi) * (1 + heta) * (1 - zeta);
                        //Ni[5, npoint] = 0.125 * (1 - ksi) * (1 + heta) * (1 - zeta);
                        //Ni[6, npoint] = 0.125 * (1 - ksi) * (1 - heta) * (1 - zeta);
                        //Ni[7, npoint] = 0.125 * (1 + ksi) * (1 - heta) * (1 - zeta);

                        Ni_ksi[0, npoint] = +0.125 * (1 + heta) * (1 + zeta); //N1_ksi(ksi,heta,zeta) sto sugkekrimeno GP "npoint" sto [,]
                        Ni_ksi[1, npoint] = -0.125 * (1 + heta) * (1 + zeta); //N2_ksi(ksi,heta,zeta) sto sugkekrimeno GP "npoint" sto [,]
                        Ni_ksi[2, npoint] = -0.125 * (1 - heta) * (1 + zeta);
                        Ni_ksi[3, npoint] = +0.125 * (1 - heta) * (1 + zeta);
                        Ni_ksi[4, npoint] = +0.125 * (1 + heta) * (1 - zeta);
                        Ni_ksi[5, npoint] = -0.125 * (1 + heta) * (1 - zeta);
                        Ni_ksi[6, npoint] = -0.125 * (1 - heta) * (1 - zeta);
                        Ni_ksi[7, npoint] = +0.125 * (1 - heta) * (1 - zeta);

                        Ni_heta[0, npoint] = 0.125 * (1 + ksi) * (+1) * (1 + zeta); //N1_heta(ksi,heta,zeta) sto sugkekrimeno GP "npoint" sto [,]
                        Ni_heta[1, npoint] = 0.125 * (1 - ksi) * (+1) * (1 + zeta); //N2_heta(ksi,heta,zeta) sto sugkekrimeno GP "npoint" sto [,]
                        Ni_heta[2, npoint] = 0.125 * (1 - ksi) * (-1) * (1 + zeta);
                        Ni_heta[3, npoint] = 0.125 * (1 + ksi) * (-1) * (1 + zeta);
                        Ni_heta[4, npoint] = 0.125 * (1 + ksi) * (+1) * (1 - zeta);
                        Ni_heta[5, npoint] = 0.125 * (1 - ksi) * (+1) * (1 - zeta);
                        Ni_heta[6, npoint] = 0.125 * (1 - ksi) * (-1) * (1 - zeta);
                        Ni_heta[7, npoint] = 0.125 * (1 + ksi) * (-1) * (1 - zeta);

                        Ni_zeta[0, npoint] = 0.125 * (1 + ksi) * (1 + heta) * (+1); //N1_zeta(ksi,heta,zeta) sto sugkekrimeno GP "npoint" sto [,]
                        Ni_zeta[1, npoint] = 0.125 * (1 - ksi) * (1 + heta) * (+1); //N2_zeta(ksi,heta,zeta) sto sugkekrimeno GP "npoint" sto [,]
                        Ni_zeta[2, npoint] = 0.125 * (1 - ksi) * (1 - heta) * (+1);
                        Ni_zeta[3, npoint] = 0.125 * (1 + ksi) * (1 - heta) * (+1);
                        Ni_zeta[4, npoint] = 0.125 * (1 + ksi) * (1 + heta) * (-1);
                        Ni_zeta[5, npoint] = 0.125 * (1 - ksi) * (1 + heta) * (-1);
                        Ni_zeta[6, npoint] = 0.125 * (1 - ksi) * (1 - heta) * (-1);
                        Ni_zeta[7, npoint] = 0.125 * (1 + ksi) * (1 - heta) * (-1);
                    }
                }
            }

            return (Ni_ksi, Ni_heta, Ni_zeta, a_123g);
        }

        // antigrafo
        public IReadOnlyList<Matrix2D> GetShapeFunctionNaturalDerivatives(IQuadrature3D quadrature)
        {
            int nGaussPoints = quadrature.IntegrationPoints.Count;
            var allNaturalDerivatives = new Matrix2D[nGaussPoints];
            for (int i = 0; i < nGaussPoints; ++i)
            {
                var gaussPoint = quadrature.IntegrationPoints[i];
                var naturalDerivatives = new double[3, 8];

                // Derivatives with respect to Xi
                naturalDerivatives[0, 0] = +0.125 * (1 + gaussPoint.Eta) * (1 + gaussPoint.Zeta);
                naturalDerivatives[0, 1] = -0.125 * (1 + gaussPoint.Eta) * (1 + gaussPoint.Zeta);
                naturalDerivatives[0, 2] = -0.125 * (1 - gaussPoint.Eta) * (1 + gaussPoint.Zeta);
                naturalDerivatives[0, 3] = +0.125 * (1 - gaussPoint.Eta) * (1 + gaussPoint.Zeta);
                naturalDerivatives[0, 4] = +0.125 * (1 + gaussPoint.Eta) * (1 - gaussPoint.Zeta);
                naturalDerivatives[0, 5] = -0.125 * (1 + gaussPoint.Eta) * (1 - gaussPoint.Zeta);
                naturalDerivatives[0, 6] = -0.125 * (1 - gaussPoint.Eta) * (1 - gaussPoint.Zeta);
                naturalDerivatives[0, 7] = +0.125 * (1 - gaussPoint.Eta) * (1 - gaussPoint.Zeta);

                // Derivatives with respect to Eta
                naturalDerivatives[1, 0] = 0.125 * (1 + gaussPoint.Xi) * (+1) * (1 + gaussPoint.Zeta);
                naturalDerivatives[1, 1] = 0.125 * (1 - gaussPoint.Xi) * (+1) * (1 + gaussPoint.Zeta);
                naturalDerivatives[1, 2] = 0.125 * (1 - gaussPoint.Xi) * (-1) * (1 + gaussPoint.Zeta);
                naturalDerivatives[1, 3] = 0.125 * (1 + gaussPoint.Xi) * (-1) * (1 + gaussPoint.Zeta);
                naturalDerivatives[1, 4] = 0.125 * (1 + gaussPoint.Xi) * (+1) * (1 - gaussPoint.Zeta);
                naturalDerivatives[1, 5] = 0.125 * (1 - gaussPoint.Xi) * (+1) * (1 - gaussPoint.Zeta);
                naturalDerivatives[1, 6] = 0.125 * (1 - gaussPoint.Xi) * (-1) * (1 - gaussPoint.Zeta);
                naturalDerivatives[1, 7] = 0.125 * (1 + gaussPoint.Xi) * (-1) * (1 - gaussPoint.Zeta);

                // Derivatives with respect to Zeta
                naturalDerivatives[2, 0] = 0.125 * (1 + gaussPoint.Xi) * (1 + gaussPoint.Eta) * (+1);
                naturalDerivatives[2, 1] = 0.125 * (1 - gaussPoint.Xi) * (1 + gaussPoint.Eta) * (+1);
                naturalDerivatives[2, 2] = 0.125 * (1 - gaussPoint.Xi) * (1 - gaussPoint.Eta) * (+1);
                naturalDerivatives[2, 3] = 0.125 * (1 + gaussPoint.Xi) * (1 - gaussPoint.Eta) * (+1);
                naturalDerivatives[2, 4] = 0.125 * (1 + gaussPoint.Xi) * (1 + gaussPoint.Eta) * (-1);
                naturalDerivatives[2, 5] = 0.125 * (1 - gaussPoint.Xi) * (1 + gaussPoint.Eta) * (-1);
                naturalDerivatives[2, 6] = 0.125 * (1 - gaussPoint.Xi) * (1 - gaussPoint.Eta) * (-1);
                naturalDerivatives[2, 7] = 0.125 * (1 + gaussPoint.Xi) * (1 - gaussPoint.Eta) * (-1);

                allNaturalDerivatives[i] = new Matrix2D(naturalDerivatives);
            }

            return allNaturalDerivatives;
        }
    }
}
