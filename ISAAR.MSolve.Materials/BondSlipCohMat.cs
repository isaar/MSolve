using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Materials.Interfaces; //using ISAAR.MSolve.PreProcessor.Interfaces;
//using ISAAR.MSolve.LinearAlgebra.Interfaces; //using ISAAR.MSolve.Matrices.Interfaces;
using ISAAR.MSolve.LinearAlgebra; //using ISAAR.MSolve.Matrices;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.LinearAlgebra.Triangulation;

namespace ISAAR.MSolve.Materials
{
    /// <summary>
    /// cohesion only friction type response in shear separation mode (with kinematic hardening like behaviour). Linear elastic behaviour in normal mode
    /// Authors Gerasimos Sotiropoulos
    /// </summary>
    public class BondSlipCohMat_v2 : ICohesiveZoneMaterial3D_v2 // TODOGerasimos
    {
        private bool modified; // opws sto MohrCoulomb gia to modified

        public double k_elastic { get; set; } // opws sto elastic 3d 
        public double k_elastic2 { get; set; }
        public double k_elastic_normal { get; set; }
        public double t_max { get; set; }
        public double[] s_0 { get; set; }
        public double[] a_0 { get; set; }
        private double[] alpha { get; set; }
        public double tol { get; set; }
        private double[] eLastConverged;
        private double[] eCurrentUpdate;
        private double[,] ConstitutiveMatrix3D;
        private double[,] ConstitutiveMatrix3Dprevious;
        private double[] stress3D;
        

        public BondSlipCohMat_v2(double k_elastic, double k_elastic2,double k_elastic_normal, double t_max, double[] s_0, double[] a_0, double tol)
        {
            this.k_elastic = k_elastic;
            this.k_elastic2 = k_elastic2;
            this.k_elastic_normal = k_elastic_normal;
            this.t_max = t_max;
            this.s_0 = s_0; //length = 2
            this.a_0 = a_0; //length = 2
            this.tol = tol;
            this.InitializeMatrices();
        }

        public BondSlipCohMat_v2(double T_o_1, double D_o_1, double k_elastic2_ratio,double T_o_3, double D_o_3, double[] s_0, double[] a_0, double tol)
        {
            this.k_elastic = T_o_1/D_o_1;
            this.k_elastic2 = k_elastic2_ratio*k_elastic;
            this.k_elastic_normal = T_o_3/D_o_3;
            this.t_max = T_o_1; // Prosoxh exei lifthei idio koino orio diarrohs fy sunolika gia th sunistamenh ths paramorfwshs anexarthtws dieftunshs
            this.s_0 = s_0; //length = 2
            this.a_0 = a_0; //length = 2
            this.tol = tol;
            this.InitializeMatrices();
        }

        ICohesiveZoneMaterial3D_v2 ICohesiveZoneMaterial3D_v2.Clone()
        {
            return this.Clone();
        }

        public BondSlipCohMat_v2 Clone()
        {
            return new BondSlipCohMat_v2(k_elastic, k_elastic2, k_elastic_normal, t_max, s_0, a_0, tol);            
        }

        private double c1;
        //private double[] sigma;

        private bool matrices_not_initialized = true;
        public void InitializeMatrices()
        {
            eCurrentUpdate = new double[2];
            eLastConverged = new double[2];// TODO: na ginetai update sto save state ennoeitai mazi me ta s_0 klp. Mporei na xrhsimopooithei h grammh apo thn arxh tou update material
            stress3D = new double[3];
            ConstitutiveMatrix3D = new double[3, 3];
            ConstitutiveMatrix3D[0, 0] = k_elastic; ConstitutiveMatrix3D[1, 1] = k_elastic; ConstitutiveMatrix3D[2, 2] = k_elastic_normal;
            ConstitutiveMatrix3Dprevious = new double[3, 3];
            ConstitutiveMatrix3Dprevious[0, 0] = k_elastic; ConstitutiveMatrix3Dprevious[1, 1] = k_elastic; ConstitutiveMatrix3Dprevious[2, 2] = k_elastic_normal;

            c1 = (k_elastic * k_elastic2) / (k_elastic - k_elastic2);
            matrices_not_initialized = false;
            //tol = Math.Pow(10, -19); //TODO: check for other convergence issues (let tol=1e-10)

            alpha = new double[2] { a_0[0], a_0[1] };
        }

        public void UpdateMaterial(double[] epsilon)
        {
            // 
            for (int k = 0; k < 3; k++)
            {
                for (int j = 0; j < 3; j++)
                {
                    ConstitutiveMatrix3Dprevious[k, j] =ConstitutiveMatrix3D[k, j];
                }
            }

            
            double[] Delta_epsilon = new double[2];
            for (int i1 = 0; i1 < 2; i1++) { Delta_epsilon[i1] = epsilon[i1] - eLastConverged[i1]; }
            eCurrentUpdate = new double[2] { epsilon[0], epsilon[1] };

            Matrix De = Matrix.CreateFromArray(new double[2, 2] { { k_elastic, 0 }, { 0, k_elastic } });
            double[] s_e = new double[2] { s_0[0] + De[0, 0] * Delta_epsilon[0], s_0[1] + De[1, 1] * Delta_epsilon[1] }; //multiplication with a diagonal matrix

            double rf = Math.Sqrt(Math.Pow((s_e[0] - a_0[0]), 2) + Math.Pow((s_e[1] - a_0[1]), 2)) - t_max;

            double Delta_l; //TODO: proswrina declare edw gia na doume ean tha kratietai
            Vector m;




            if (rf < 0)
            {
                stress3D = new double[3] { s_e[0], s_e[1], k_elastic_normal * epsilon[2] }; //oi duo prwtoi oroi einai to "sigma"
                alpha = new double[2] { a_0[0], a_0[1] };
                Delta_l = 0;
                m =  Vector.CreateZero(2);
                ConstitutiveMatrix3D = new double[3, 3];
                ConstitutiveMatrix3D[0, 0] = k_elastic; ConstitutiveMatrix3D[1, 1] = k_elastic; ConstitutiveMatrix3D[2, 2] = k_elastic_normal;                
            }
            else
            {
                Delta_l = 0; // arxikh timh 
                double[] sigma = new double[2] { s_e[0], s_e[1] };
                alpha = new double[2] { a_0[0], a_0[1] };

                double[] res_vec = new double[5];
                res_vec[4] = rf;

                double com_term1 = Math.Sqrt(Math.Pow(sigma[0] - alpha[0], 2) + Math.Pow(sigma[1] - alpha[1], 2));

                m = Vector.CreateZero(2);
                for (int i1 = 0; i1 < 2; i1++) { m[i1] = (1 / com_term1) * (sigma[i1] - alpha[i1]); }

                Matrix tangent_matrix =Matrix.CreateZero(5,5);
                Matrix tangent_matrix_inv=Matrix.CreateZero(5, 5);
                int iter = 0;
                while (((Vector.CreateFromArray(res_vec)).Norm2() / t_max) > tol || iter == 0)
                {
                    //compute tangential matrix (residual derivatives)
                    double com_term2 = Math.Pow(Math.Pow(sigma[0] - alpha[0], 2) + Math.Pow(sigma[1] - alpha[1], 2), (-((double)3 /(double)2)));
                    Matrix dmi_dti = Matrix.CreateZero(2, 2);
                    dmi_dti[0, 0] = com_term2 * Math.Pow(sigma[1] - alpha[1], 2);
                    dmi_dti[0, 1] = -com_term2 * (sigma[0] - alpha[0]) * (sigma[1] - alpha[1]);
                    dmi_dti[1, 0] = dmi_dti[0, 1];
                    dmi_dti[1, 1] = com_term2 * Math.Pow(sigma[0] - alpha[0], 2);

                    Matrix dmi_dai = Matrix.CreateZero(2, 2);
                    dmi_dai[0, 0] = -dmi_dti[0, 0];
                    dmi_dai[0, 1] = -dmi_dti[0, 1];
                    dmi_dai[1, 0] = -dmi_dti[1, 0];
                    dmi_dai[1, 1] = -dmi_dti[1, 1];

                    Matrix drs_ds = (De * dmi_dti);
                    drs_ds.ScaleIntoThis(Delta_l);
                    drs_ds[0, 0] += 1;
                    drs_ds[1, 1] += 1;

                    Matrix drs_da = (De * dmi_dai);
                    drs_da.ScaleIntoThis(Delta_l);

                    double[] drs_dl = new double[2];
                    for (int i1 = 0; i1 < 2; i1++) { for (int i2 = 0; i2 < 2; i2++) { drs_dl[i1] += De[i1, i2] * m[i2]; } }


                    Matrix dra_ds = Matrix.CreateFromArray(new double[2, 2] { { dmi_dti[0, 0], dmi_dti[0, 1] }, { dmi_dti[1, 0], dmi_dti[1, 1] } });
                    dra_ds.ScaleIntoThis(-Delta_l * c1);

                    Matrix dra_dra =Matrix.CreateFromArray(new double[2, 2] { { dmi_dai[0, 0], dmi_dai[0, 1] }, { dmi_dai[1, 0], dmi_dai[1, 1] } });
                    dra_dra.ScaleIntoThis(-Delta_l * c1);
                    dra_dra[0, 0] += 1;
                    dra_dra[1, 1] += 1; //TODO check new structures matrix for this

                    double[] dra_dl = new double[2] { -c1 * m[0], -c1 * m[1] };

                    double[] drf_ds = new double[2] { m[0], m[1] };

                    double drf_da1 = (1 / com_term1) * (-sigma[0] + alpha[0]);
                    double drf_da2 = (1 / com_term1) * (-sigma[1] + alpha[1]);

                    double drf_dl = 0;

                    tangent_matrix = Matrix.CreateFromArray(new double[5, 5]{ {drs_ds[0,0],drs_ds[0,1],drs_da[0,0],drs_da[0,1],drs_dl[0]},
                        {drs_ds[1,0],drs_ds[1,1],drs_da[1,0],drs_da[1,1],drs_dl[1]},
                        {dra_ds[0,0],dra_ds[0,1],dra_dra[0,0],dra_dra[0,1],dra_dl[0]},
                        {dra_ds[1,0],dra_ds[1,1],dra_dra[1,0],dra_dra[1,1],dra_dl[1]},
                        {drf_ds[0],drf_ds[1],drf_da1,drf_da2,drf_dl}});

                    // solution and update of vectors
                    tangent_matrix_inv = tangent_matrix.Invert();
                    Vector Solution = tangent_matrix_inv * (Vector.CreateFromArray(res_vec));//Vector Solution = tangent_matrix.SolveLU(res_vec, false);
                    sigma[0] += -Solution[0];
                    sigma[1] += -Solution[1];
                    alpha[0] += -Solution[2];
                    alpha[1] += -Solution[3];
                    Delta_l += -Solution[4];

                    com_term1 = Math.Sqrt(Math.Pow(sigma[0] - alpha[0], 2) + Math.Pow(sigma[1] - alpha[1], 2));

                    for (int i2 = 0; i2 < 2; i2++) { m[i2] = (1 / com_term1) * (sigma[i2] - alpha[i2]); }

                    double[] rs = new double[2];
                    for (int i1 = 0; i1 < 2; i1++) { for (int i2 = 0; i2 < 2; i2++) { rs[i1] += Delta_l * De[i1, i2] * m[i2]; } }
                    for (int i2 = 0; i2 < 2; i2++) { rs[i2] += sigma[i2] - s_e[i2]; }

                    double[] ra = new double[2];
                    for (int i2 = 0; i2 < 2; i2++) { ra[i2] += alpha[i2] - a_0[i2] - Delta_l * c1 * m[i2]; }

                    rf = Math.Sqrt(Math.Pow(sigma[0] - alpha[0], 2) + Math.Pow(sigma[1] - alpha[1], 2)) - t_max;

                    res_vec = new double[5] { rs[0], rs[1], ra[0], ra[1], rf };
                    iter += 1;
                }                
                //var tangent_matrix_inv = tangent_matrix.SolveLU(eye5, false);

                for (int i1 = 0; i1 < 2; i1++)
                {
                    for (int i2 = 0; i2 < 2; i2++)
                    {
                        ConstitutiveMatrix3D[i1, i2] = 0;
                        for (int i3 = 0; i3 < 2; i3++) { ConstitutiveMatrix3D[i1, i2] += tangent_matrix_inv[i1, i3] * De[i3, i2]; }
                    }
                }

                stress3D = new double[3] { sigma[0], sigma[1], k_elastic_normal * epsilon[2] };


            }

            this.modified = CheckIfConstitutiveMatrixChanged();
        }

        private bool CheckIfConstitutiveMatrixChanged()
        {
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    if (Math.Abs(ConstitutiveMatrix3Dprevious[i, j] - ConstitutiveMatrix3D[i, j]) > 1e-10)
                        return true;

            return false;
        }

        public double[] Tractions // opws xrhsimopoeitai sto mohrcoulomb kai hexa8
        {
            get { return stress3D; }
        }

        public IMatrixView ConstitutiveMatrix
        {
            get
            {

                return  Matrix.CreateFromArray(ConstitutiveMatrix3D);
            }
        }

        public void SaveState()
        {
            a_0 = new double[2] { alpha[0], alpha[1] };
            s_0 = new double[2] { stress3D[0], stress3D[1] };
            eLastConverged = new double[2] { eCurrentUpdate[0], eCurrentUpdate[1] };
        }

        public bool Modified
        {
            get { return modified; }
        }

        public void ResetModified()
        {
            modified = false;
        }

        public int ID
        {
            get { return 1000; }
        }

        public void ClearState() // pithanws TODO 
        {
            //ean thelei to D_tan ths arxikhs katastashs tha epistrepsoume const me De
            // alla oxi ia na to xrhsimopoihsei gia elastiko se alles periptwseis
            //opws
            // sthn epanalhptikh diadikasia (opws px provider.Reset pou sumvainei se polles epanalipseis?)
        }
        public void ClearTractions()
        {

        }

        private double youngModulus = 1;
        public double YoungModulus
        {
            get { throw new InvalidOperationException(); }
            set { throw new InvalidOperationException(); }
        }

        private double poissonRatio = 1;
        public double PoissonRatio
        {
            get { return poissonRatio; }
            set { throw new InvalidOperationException(); }
        }

        private double[] coordinates;
        public double[] Coordinates
        {

            get { return coordinates; }
            set { throw new InvalidOperationException(); }
        }
    }
}
