using System;
using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.Materials.Interfaces;

namespace ISAAR.MSolve.FEM.Elements
{
    public class Hexa8Memoizer
    {
        private readonly Dictionary<int, Tuple<double[], double[, ,]>> integrationDictionary = new Dictionary<int, Tuple<double[], double[, ,]>>();

        public Tuple<double[], double[, ,]> GetIntegrationData(int element)
        {
            if (integrationDictionary.ContainsKey(element))
                return integrationDictionary[element];
            else
                return new Tuple<double[],double[,,]>(null, null);
        }

        public void SetIntegrationData(int element, Tuple<double[], double[, ,]> integrationData)
        {
            integrationDictionary.Add(element, integrationData);
        }
    }

    public class Hexa8WithStochasticMaterial : Hexa8
    {
        protected readonly new IStochasticContinuumMaterial3D[] materialsAtGaussPoints;
        protected readonly Hexa8Memoizer memoizer;

        private static double[][] integrationPoints = new double[][] 
        { 
            new double[] { }, 
            new double[] { 0 }, 
            new double[] { -0.5773502691896, 0.5773502691896 }, 
            new double[] { -0.774596669, 0, 0.774596669 }, 
            new double[] { -0.8611363115941, -0.3399810435849, 0.3399810435849, 0.8611363115941 } 
        };

        public Hexa8WithStochasticMaterial(IStochasticContinuumMaterial3D material)
        {
            //materialsAtGaussPoints = new IStochasticContinuumMaterial3D[iInt3];
            //for (int i = 0; i < iInt3; i++)
            //    materialsAtGaussPoints[i] = (IStochasticContinuumMaterial3D)material.Clone();
        }

        public Hexa8WithStochasticMaterial(IStochasticContinuumMaterial3D material, Hexa8Memoizer memoizer) : this(material)
        {
            this.memoizer = memoizer;
        }

        //public Hexa8WithStochasticMaterial(IFiniteElementMaterial3D material, IStochasticCoefficientsProvider coefficientsProvider)
        //    : this(material)
        //{
        //    this.coefficientsProvider = coefficientsProvider;
        //}

        public override IMatrix2D StiffnessMatrix(IElement element)
        {
            double[, ,] afE = new double[iInt3, 6, 6];
            int iPos = 0;
            for (int i1 = 0; i1 < iInt; i1++)
                for (int i2 = 0; i2 < iInt; i2++)
                    for (int i3 = 0; i3 < iInt; i3++)
                    {
                        iPos = i1 * iInt2 + i2 * iInt + i3;
                        var e = ((Matrix2D)materialsAtGaussPoints[iPos].GetConstitutiveMatrix(GetStochasticPoints(element, i1, i2, i3)));
                        for (int j = 0; j < 6; j++)
                            for (int k = 0; k < 6; k++)
                                afE[iPos, j, k] = e[j, k];
                        //afE[i, j, k] = ((Matrix2D<double>)materialsAtGaussPoints[i].GetConstitutiveMatrix(GetStochasticPoints(element, i / iInt2, (i % iInt2) / iInt, i % iInt)))[j, k];
                    }

            double[, ,] faB;
            double[] faWeight;
            Tuple<double[], double[, ,]> integrationData = new Tuple<double[],double[,,]>(null, null);
            if (memoizer != null)
                integrationData = memoizer.GetIntegrationData(element.ID);
            if (integrationData.Item1 == null)
            {
                faB = new double[iInt3, 24, 6];
                faWeight = new double[iInt3];
                double[,] faXYZ = GetCoordinates(element);
                double[,] faDS = new double[iInt3, 24];
                double[,] faS = new double[iInt3, 8];
                double[] faDetJ = new double[iInt3];
                double[, ,] faJ = new double[iInt3, 3, 3];
                CalcH8GaussMatrices(ref iInt, faXYZ, faWeight, faS, faDS, faJ, faDetJ, faB);
                //for (int i = 0; i < iInt; i++)
                //    for (int j = 0; j < iInt; j++)
                //        for (int k = 0; k < iInt; k++)
                //        {
                //            faWeight[i * iInt * iInt + j * iInt + k] *= coefficientsProvider.GetCoefficient(GetStochasticPoints(element, i, j, k));
                //            //faWeight[i * iInt * iInt + j * iInt + k] *= coefficientsProvider.GetCoefficient(new double[] { integrationPoints[iInt][i], integrationPoints[iInt][j], integrationPoints[iInt][k] });
                //        }
                if (memoizer != null)
                    memoizer.SetIntegrationData(element.ID, new Tuple<double[], double[, ,]>(faWeight, faB));
            }
            else
            {
                faB = integrationData.Item2;
                faWeight = integrationData.Item1;
            }
            double[] faK = new double[300];
            CalcH8K(ref iInt, afE, faB, faWeight, faK);
            return dofEnumerator.GetTransformedMatrix(new SymmetricMatrix2D(faK));
        }

        private double[] GetStochasticPoints(IElement element, int iX, int iY, int iZ)
        {
            // Calculate for element centroid
            double X = 0;
            double Y = 0;
            double Z = 0;
            double minX = element.INodes[0].X;
            double minY = element.INodes[0].Y;
            double minZ = element.INodes[0].Z;

            for (int i = 0; i < 8; i++)
            {
                minX = minX > element.INodes[i].X ? element.INodes[i].X : minX;
                minY = minY > element.INodes[i].Y ? element.INodes[i].Y : minY;
                minZ = minZ > element.INodes[i].Z ? element.INodes[i].Z : minZ;
                for (int j = i + 1; j < 8; j++)
                {
                    X = X < Math.Abs(element.INodes[j].X - element.INodes[i].X) ? Math.Abs(element.INodes[j].X - element.INodes[i].X) : X;
                    Y = Y < Math.Abs(element.INodes[j].Y - element.INodes[i].Y) ? Math.Abs(element.INodes[j].Y - element.INodes[i].Y) : Y;
                    Z = Z < Math.Abs(element.INodes[j].Z - element.INodes[i].Z) ? Math.Abs(element.INodes[j].Z - element.INodes[i].Z) : Z;
                }
            }

            double pointX = minX + X / 2;
            double pointY = minY + Y / 2;
            double pointZ = minZ + Z / 2;

            return new double[] { pointX, pointY, pointZ };

            //// Calculate for individual gauss point
            ////if (iInt != 2) throw new ArgumentException("Stochastic provided functions with integration order of 2.");

            //double X = 0;
            //double Y = 0;
            //double Z = 0;
            //double minX = element.Nodes[0].X;
            //double minY = element.Nodes[0].Y;
            //double minZ = element.Nodes[0].Z;

            //for (int i = 0; i < 8; i++)
            //{
            //    minX = minX > element.Nodes[i].X ? element.Nodes[i].X : minX;
            //    minY = minY > element.Nodes[i].Y ? element.Nodes[i].Y : minY;
            //    minZ = minZ > element.Nodes[i].Z ? element.Nodes[i].Z : minZ;
            //    for (int j = i + 1; j < 8; j++)
            //    {
            //        X = X < Math.Abs(element.Nodes[j].X - element.Nodes[i].X) ? Math.Abs(element.Nodes[j].X - element.Nodes[i].X) : X;
            //        Y = Y < Math.Abs(element.Nodes[j].Y - element.Nodes[i].Y) ? Math.Abs(element.Nodes[j].Y - element.Nodes[i].Y) : Y;
            //        Z = Z < Math.Abs(element.Nodes[j].Z - element.Nodes[i].Z) ? Math.Abs(element.Nodes[j].Z - element.Nodes[i].Z) : Z;
            //    }
            //}

            //double pointX = minX + X / 2 * (integrationPoints[iInt][iX] + 1);
            //double pointY = minY + Y / 2 * (integrationPoints[iInt][iY] + 1);
            //double pointZ = minZ + Z / 2 * (integrationPoints[iInt][iZ] + 1);

            //return new double[] { pointX, pointY, pointZ };
        }
    }
}
