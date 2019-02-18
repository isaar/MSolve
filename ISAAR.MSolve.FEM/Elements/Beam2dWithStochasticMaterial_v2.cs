using System;
using System.Collections.Generic;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Materials;
using ISAAR.MSolve.LinearAlgebra.Matrices;
using ISAAR.MSolve.Materials.Interfaces;

namespace ISAAR.MSolve.FEM.Elements
{
    public class Beam2DMemoizer
    {
        private readonly Dictionary<int, Tuple<double[], double[,,]>> integrationDictionary =
            new Dictionary<int, Tuple<double[], double[,,]>>();

        public Tuple<double[], double[,,]> GetIntegrationData(int element)
        {
            if (integrationDictionary.ContainsKey(element))
                return integrationDictionary[element];
            else
                return new Tuple<double[], double[,,]>(null, null);
        }

        public void SetIntegrationData(int element, Tuple<double[], double[,,]> integrationData)
        {
            integrationDictionary.Add(element, integrationData);
        }
    }

    public class Beam2DWithStochasticMaterial_v2 : EulerBeam2D_v2
    {
        protected readonly new IStochasticContinuumMaterial3D_v2 Material;
        protected readonly Beam2DMemoizer memoizer;

        public Beam2DWithStochasticMaterial_v2(IStochasticContinuumMaterial3D_v2 material): 
            base((material as StochasticElasticMaterial).YoungModulus)
        {
            this.Material = material;
        }

        public Beam2DWithStochasticMaterial_v2(IStochasticContinuumMaterial3D_v2 material, Beam2DMemoizer memoizer) :
            this(material)
        {
            this.memoizer = memoizer;
        }

        //public Hexa8WithStochasticMaterial(IFiniteElementMaterial3D material, IStochasticCoefficientsProvider coefficientsProvider)
        //    : this(material)
        //{
        //    this.coefficientsProvider = coefficientsProvider;
        //}

        public override IMatrix StiffnessMatrix(IElement_v2 element)
        {
            double x2 = Math.Pow(element.Nodes[1].X - element.Nodes[0].X, 2);
            double y2 = Math.Pow(element.Nodes[1].Y - element.Nodes[0].Y, 2);
            double L = Math.Sqrt(x2 + y2);
            double c = (element.Nodes[1].X - element.Nodes[0].X) / L;
            double c2 = c * c;
            double s = (element.Nodes[1].Y - element.Nodes[0].Y) / L;
            double s2 = s * s;
            double[] coordinates = GetStochasticPoints(element);
            double EL = (Material as StochasticElasticMaterial).GetStochasticMaterialProperties(coordinates)[0] / L;
            double EAL = EL * SectionArea;
            double EIL = EL * MomentOfInertia;
            double EIL2 = EIL / L;
            double EIL3 = EIL2 / L;
            int order = 6;
            return DofEnumerator.GetTransformedMatrix(SymmetricMatrix.CreateFromPackedRowMajorArray(new double[]
            {
                c2 * EAL + 12 * s2 * EIL3, c * s * EAL - 12 * c * s * EIL3, -6 * s * EIL2, -c2 * EAL - 12 * s2 * EIL3, -c * s * EAL + 12 * c * s * EIL3, -6 * s * EIL2,
                s2 * EAL + 12 * c2 * EIL3, 6 * c * EIL2, -s * c * EAL + 12 * c * s * EIL3, -s2 * EAL - 12 * c2 * EIL3,  6 * c * EIL2,
                4 * EIL, 6 * s * EIL2, -6 * c * EIL2, 2 * EIL,
                c2 * EAL + 12 * s2 * EIL3, s * c * EAL - 12 * c * s * EIL3, 6 * s * EIL2,
                s2 * EAL + 12 * c2 * EIL3, -6 * c * EIL2,
                4 * EIL
            }, order));
        }

        private double[] GetStochasticPoints(IElement_v2 element)
        {
            // Calculate for element centroid
            double X = 0;
            double Y = 0;
            double minX = element.Nodes[0].X;
            double minY = element.Nodes[0].Y;

            for (int i = 0; i < 2; i++)
            {
                minX = minX > element.Nodes[i].X ? element.Nodes[i].X : minX;
                minY = minY > element.Nodes[i].Y ? element.Nodes[i].Y : minY;
                for (int j = i + 1; j < 2; j++)
                {
                    X = X < Math.Abs(element.Nodes[j].X - element.Nodes[i].X)
                        ? Math.Abs(element.Nodes[j].X - element.Nodes[i].X)
                        : X;
                    Y = Y < Math.Abs(element.Nodes[j].Y - element.Nodes[i].Y)
                        ? Math.Abs(element.Nodes[j].Y - element.Nodes[i].Y)
                        : Y;
                }
            }

            double pointX = minX + X / 2;
            double pointY = minY + Y / 2;

            return new double[] {pointX, pointY};

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