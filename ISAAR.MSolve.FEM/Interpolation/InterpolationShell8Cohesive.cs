using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.FEM.Interpolation
{
    class InterpolationShell8Cohesive
    {
        private static readonly InterpolationShell8Cohesive uniqueInstance = new InterpolationShell8Cohesive();

        public static InterpolationShell8Cohesive UniqueInstance => uniqueInstance;

        public (double[][] N1, double[][,] N3, double[][] N1_ksi, double[][] N1_heta, double[] a_12g) GetShapeFunctionAndGaussPointData(int nGaussPoints, int gp_d1_coh ,int gp_d2_coh)
        {
            //initialize edw twn arrays pou pleon tha epistrefontai
            double[] a_12g;
            double[][] N1;
            double[][,] N3;
            double[][] N1_ksi;
            double[][] N1_heta;

            // PROSTHIKI RAM
            double a_1g = 0;
            double a_2g = 0;
            double[] N_i;
            double[] N_i_ksi;
            double[] N_i_heta;
            double ksi = 0;
            double heta = 0;
            //double zeta;
            int npoint;


            nGaussPoints = gp_d1_coh * gp_d2_coh;
            a_12g = new double[nGaussPoints];
            N_i = new double[8]; // 4 gia ligoterous komvous coh8
            N_i_ksi = new double[8]; // 
            N_i_heta = new double[8]; // 
            N1 = new double[nGaussPoints][];
            N3 = new double[nGaussPoints][,];
            N1_ksi = new double[nGaussPoints][];
            N1_heta = new double[nGaussPoints][];
            for (int l = 0; l < nGaussPoints; l++)
            {
                N1[l] = new double[8];
                N3[l] = new double[3, 24];
                N1_ksi[l] = new double[8];
                N1_heta[l] = new double[8];
            }
            for (int j = 0; j < gp_d1_coh; j++)
            {
                for (int k = 0; k < gp_d2_coh; k++)
                {
                    npoint = j * gp_d1_coh + k;
                    if (gp_d1_coh == 3)
                    {
                        ksi = 0.5 * (j - 1) * (j - 2) * (-0.774596669241483) + (-1) * (j) * (j - 2) * (0) + 0.5 * (j) * (j - 1) * (0.774596669241483);
                        a_1g = 0.5 * (j - 1) * (j - 2) * (0.555555555555555) + (-1) * (j) * (j - 2) * (0.888888888888888) + 0.5 * (j) * (j - 1) * (0.555555555555555);
                    }
                    if (gp_d1_coh == 2)
                    {
                        ksi = (-0.577350269189626) * (j - 1) * (-1) + (0.577350269189626) * (j) * (+1);
                        a_1g = 1;
                    }
                    if (gp_d2_coh == 3)
                    {
                        heta = 0.5 * (k - 1) * (k - 2) * (-0.774596669241483) + (-1) * (k) * (k - 2) * (0) + 0.5 * (k) * (k - 1) * (0.774596669241483);
                        a_2g = 0.5 * (k - 1) * (k - 2) * (0.555555555555555) + (-1) * (k) * (k - 2) * (0.888888888888888) + 0.5 * (k) * (k - 1) * (0.555555555555555);
                    }
                    if (gp_d2_coh == 2)
                    {
                        heta = (-0.577350269189626) * (k - 1) * (-1) + (0.577350269189626) * (k) * (+1);
                        a_2g = 1;
                    }

                    a_12g[npoint] = a_1g * a_2g;

                    N_i[4] = 0.5 * (1 - Math.Pow(ksi, 2)) * (1 + heta);
                    N_i[5] = 0.5 * (1 - Math.Pow(heta, 2)) * (1 - ksi);
                    N_i[6] = 0.5 * (1 - Math.Pow(ksi, 2)) * (1 - heta);
                    N_i[7] = 0.5 * (1 - Math.Pow(heta, 2)) * (1 + ksi);

                    N_i[0] = 0.25 * (1 + ksi) * (1 + heta) - 0.5 * N_i[4] - 0.5 * N_i[7];
                    N_i[1] = 0.25 * (1 - ksi) * (1 + heta) - 0.5 * N_i[4] - 0.5 * N_i[5];
                    N_i[2] = 0.25 * (1 - ksi) * (1 - heta) - 0.5 * N_i[5] - 0.5 * N_i[6];
                    N_i[3] = 0.25 * (1 + ksi) * (1 - heta) - 0.5 * N_i[6] - 0.5 * N_i[7];

                    N_i_ksi[4] = (-ksi) * (1 + heta);
                    N_i_ksi[5] = -0.5 * (1 - Math.Pow(heta, 2));
                    N_i_ksi[6] = 0.5 * (-2 * ksi) * (1 - heta);
                    N_i_ksi[7] = 0.5 * (1 - Math.Pow(heta, 2));
                    N_i_ksi[0] = +0.25 * (1 + heta) - 0.5 * N_i_ksi[4] - 0.5 * N_i_ksi[7];
                    N_i_ksi[1] = -0.25 * (1 + heta) - 0.5 * N_i_ksi[4] - 0.5 * N_i_ksi[5];
                    N_i_ksi[2] = -0.25 * (1 - heta) - 0.5 * N_i_ksi[5] - 0.5 * N_i_ksi[6];
                    N_i_ksi[3] = +0.25 * (1 - heta) - 0.5 * N_i_ksi[6] - 0.5 * N_i_ksi[7];

                    N_i_heta[4] = 0.5 * (1 - Math.Pow(ksi, 2));
                    N_i_heta[5] = 0.5 * (-2 * heta) * (1 - ksi);
                    N_i_heta[6] = 0.5 * (1 - Math.Pow(ksi, 2)) * (-1);
                    N_i_heta[7] = 0.5 * (-2 * heta) * (1 + ksi);
                    N_i_heta[0] = +0.25 * (1 + ksi) - 0.5 * N_i_heta[4] - 0.5 * N_i_heta[7];
                    N_i_heta[1] = +0.25 * (1 - ksi) - 0.5 * N_i_heta[4] - 0.5 * N_i_heta[5];
                    N_i_heta[2] = -0.25 * (1 - ksi) - 0.5 * N_i_heta[5] - 0.5 * N_i_heta[6];
                    N_i_heta[3] = -0.25 * (1 + ksi) - 0.5 * N_i_heta[6] - 0.5 * N_i_heta[7];

                    for (int l = 0; l < 8; l++)  // to 8 ginetai 4 gia to cohesive8node
                    { N1[npoint][l] = N_i[l]; }

                    for (int l = 0; l < 3; l++)  // arxika mhdenismos twn stoixweiwn tou pinaka
                    {
                        for (int m = 0; m < 24; m++)
                        { N3[npoint][l, m] = 0; }
                    }

                    for (int l = 0; l < 3; l++)
                    {
                        for (int m = 0; m < 8; m++)
                        { N3[npoint][l, l + 3 * m] = N_i[m]; }
                    }

                    for (int l = 0; l < 8; l++)
                    { N1_ksi[npoint][l] = N_i_ksi[l]; }

                    for (int l = 0; l < 8; l++)
                    { N1_heta[npoint][l] = N_i_heta[l]; }

                }
            }

            return (N1, N3, N1_ksi, N1_heta, a_12g);
            
        }

    }
}
