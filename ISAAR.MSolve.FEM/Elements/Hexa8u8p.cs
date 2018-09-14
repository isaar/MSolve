using System;
using System.Collections.Generic;
using System.Runtime.InteropServices;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.Materials.Interfaces;

namespace ISAAR.MSolve.FEM.Elements
{
    public class Hexa8u8p : IPorousFiniteElement
    {
        protected static int iInt = 2;
        protected static int iInt2 = iInt * iInt;
        protected static int iInt3 = iInt * iInt * iInt;
        private readonly static DOFType[] nodalDOFTypes = new DOFType[] { DOFType.X, DOFType.Y, DOFType.Z, DOFType.Pore };
        private readonly static DOFType[][] dofTypes = new DOFType[][] { nodalDOFTypes, nodalDOFTypes, nodalDOFTypes,
            nodalDOFTypes, nodalDOFTypes, nodalDOFTypes, nodalDOFTypes, nodalDOFTypes };
        protected IIsotropicContinuumMaterial3D[] materialsAtGaussPoints;
        protected IElementDOFEnumerator dofEnumerator = new GenericDOFEnumerator();

        #region Fortran imports
        [DllImport("femelements.dll",
            EntryPoint = "CALCH8GAUSSMATRICES",
            CallingConvention = CallingConvention.Cdecl)]
        protected static extern void CalcH8GaussMatrices(ref int iInt, [MarshalAs(UnmanagedType.LPArray)]double[,] faXYZ,
            [MarshalAs(UnmanagedType.LPArray)]double[] faWeight, [MarshalAs(UnmanagedType.LPArray)]double[,] faS,
            [MarshalAs(UnmanagedType.LPArray)]double[,] faDS, [MarshalAs(UnmanagedType.LPArray)]double[, ,] faJ,
            [MarshalAs(UnmanagedType.LPArray)]double[] faDetJ, [MarshalAs(UnmanagedType.LPArray)]double[, ,] faB);

        [DllImport("femelements.dll",
            EntryPoint = "CALCH8STRAINS",
            CallingConvention = CallingConvention.Cdecl)]
        private static extern void CalcH8Strains(ref int iInt,
            [MarshalAs(UnmanagedType.LPArray)]double[, ,] faB, [MarshalAs(UnmanagedType.LPArray)]double[] fau,
            [MarshalAs(UnmanagedType.LPArray)]double[,] faStrains);

        [DllImport("femelements.dll",
            EntryPoint = "CALCH8FORCES",
            CallingConvention = CallingConvention.Cdecl)]
        private static extern void CalcH8Forces(ref int iInt,
            [MarshalAs(UnmanagedType.LPArray)]double[, ,] faB, [MarshalAs(UnmanagedType.LPArray)]double[] faWeight,
            [MarshalAs(UnmanagedType.LPArray)]double[,] faStresses,
            [MarshalAs(UnmanagedType.LPArray)]double[] faForces);

        [DllImport("femelements.dll",
            EntryPoint = "CALCH20U8PFORCESWATERACC",
            CallingConvention = CallingConvention.Cdecl)]
        private static extern void CalcH20u8pForcesWaterAcc(ref int iInt,
            [MarshalAs(UnmanagedType.LPArray)]bool[] alImpermeable, ref double ffDensity, 
            [MarshalAs(UnmanagedType.LPArray)]double[] afPermeability,
            [MarshalAs(UnmanagedType.LPArray)]double[] afXw,
            [MarshalAs(UnmanagedType.LPArray)]double[] afSw,
            [MarshalAs(UnmanagedType.LPArray)]double[] afAcc,
            [MarshalAs(UnmanagedType.LPArray)]double[,] afS,
            [MarshalAs(UnmanagedType.LPArray)]double[, ,] afB, [MarshalAs(UnmanagedType.LPArray)]double[] afWeights,
            [MarshalAs(UnmanagedType.LPArray)]double[] afLocalForces);

        [DllImport("femelements.dll",
            EntryPoint = "CALCH8K",
            CallingConvention = CallingConvention.Cdecl)]
        protected static extern void CalcH8K(ref int iInt, [MarshalAs(UnmanagedType.LPArray)]double[, ,] faE,
            [MarshalAs(UnmanagedType.LPArray)]double[, ,] faB, [MarshalAs(UnmanagedType.LPArray)]double[] faWeight,
            [MarshalAs(UnmanagedType.LPArray)]double[] faK);

        [DllImport("femelements.dll",
            EntryPoint = "CALCH8MLUMPED",
            CallingConvention = CallingConvention.Cdecl)]
        private static extern void CalcH8MLumped(ref int iInt, ref double fDensity,
            [MarshalAs(UnmanagedType.LPArray)]double[] faWeight, [MarshalAs(UnmanagedType.LPArray)]double[] faM);

        [DllImport("femelements.dll",
            EntryPoint = "CALCH20U8PH",
            CallingConvention = CallingConvention.Cdecl)]
        private static extern void CalcH20u8pH(ref int iInt, [MarshalAs(UnmanagedType.LPArray)]double[] faPermeability,
            [MarshalAs(UnmanagedType.LPArray)]double[,] faS, [MarshalAs(UnmanagedType.LPArray)]double[, ,] faB,
            [MarshalAs(UnmanagedType.LPArray)]double[] faWeight, [MarshalAs(UnmanagedType.LPArray)]double[] faH);

        [DllImport("femelements.dll",
            EntryPoint = "CALCH20U8PS",
            CallingConvention = CallingConvention.Cdecl)]
        private static extern void CalcH20u8pS(ref int iInt, [MarshalAs(UnmanagedType.LPArray)]double[] faXwDivQ,
            [MarshalAs(UnmanagedType.LPArray)]double[,] faS,
            [MarshalAs(UnmanagedType.LPArray)]double[] faWeight, [MarshalAs(UnmanagedType.LPArray)]double[] faM);

        [DllImport("femelements.dll",
            EntryPoint = "CALCH8U8PQMINUS",
            CallingConvention = CallingConvention.Cdecl)]
        private static extern void CalcH8u8pQMinus(ref int iInt, ref double fPoreA, 
            [MarshalAs(UnmanagedType.LPArray)]double[] faXw,
            [MarshalAs(UnmanagedType.LPArray)]double[, ,] faB, [MarshalAs(UnmanagedType.LPArray)]double[,] faS,
            [MarshalAs(UnmanagedType.LPArray)]double[] faWeight, [MarshalAs(UnmanagedType.LPArray)]double[,] faQ);

        #endregion

        protected Hexa8u8p()
        {
        }

        public Hexa8u8p(IIsotropicContinuumMaterial3D material)
        {
            //materialsAtGaussPoints = new IIsotropicContinuumMaterial3D[Hexa8u8p.iInt3];
            //for (int i = 0; i < Hexa8u8p.iInt3; i++)
            //    materialsAtGaussPoints[i] = (IIsotropicContinuumMaterial3D)material.Clone();
        }

        public Hexa8u8p(IIsotropicContinuumMaterial3D material, IElementDOFEnumerator dofEnumerator) : this(material)
        {
            this.dofEnumerator = dofEnumerator;
        }

        public IElementDOFEnumerator DOFEnumerator
        {
            get { return dofEnumerator; }
            set { dofEnumerator = value; }
        }

        public double FluidBulkModulus { get; set; }
        public double SolidDensity { get; set; }
        public double FluidDensity { get; set; }
        public double Permeability { get; set; }
        public double Porosity { get; set; }
        public double Saturation { get; set; }
        public double PoreA { get; set; }
        public double Xw { get; set; }
        public double Cs { get; set; }
        public double RayleighAlpha { get; set; }
        public double RayleighBeta { get; set; }

        public double SolidBulkModulus { get { return materialsAtGaussPoints[0].YoungModulus / (3 - 6 * materialsAtGaussPoints[0].PoissonRatio); } }
        public double Density { get { return Porosity * Saturation * FluidDensity + 
            (1 - Porosity) * SolidDensity; } }
        public double QInv { get { return Cs + Porosity * Saturation / FluidBulkModulus + 
            (PoreA - Porosity) * Saturation / SolidBulkModulus; } }

        protected double[,] GetCoordinates(IElement element)
        {
            double[,] faXYZ = new double[dofTypes.Length, 3];
            for (int i = 0; i < dofTypes.Length; i++)
            {
                faXYZ[i, 0] = element.INodes[i].X;
                faXYZ[i, 1] = element.INodes[i].Y;
                faXYZ[i, 2] = element.INodes[i].Z;
            }
            return faXYZ;
        }

        #region IElementType Members

        public int ID
        {
            get { return 11; }
        }

        public ElementDimensions ElementDimensions
        {
            get { return ElementDimensions.ThreeD; }
        }

        public IList<IList<DOFType>> GetElementDOFTypes(IElement element)
        {
            return dofTypes;
        }

        public IList<Node> GetNodesForMatrixAssembly(Element element)
        {
            return element.Nodes;
        }

        public virtual IMatrix2D StiffnessMatrix(IElement element)
        {
            double[, ,] afE = new double[iInt3, 6, 6];
            for (int i = 0; i < iInt3; i++)
                for (int j = 0; j < 6; j++)
                    for (int k = 0; k < 6; k++)
                        afE[i, j, k] = ((Matrix2D)materialsAtGaussPoints[i].ConstitutiveMatrix)[j, k];
            double[,] faXYZ = GetCoordinates(element);
            double[,] faDS = new double[iInt3, 24];
            double[,] faS = new double[iInt3, 8];
            double[, ,] faB = new double[iInt3, 24, 6];
            double[] faDetJ = new double[iInt3];
            double[, ,] faJ = new double[iInt3, 3, 3];
            double[] faWeight = new double[iInt3];
            double[] faK = new double[300];
            CalcH8GaussMatrices(ref iInt, faXYZ, faWeight, faS, faDS, faJ, faDetJ, faB);
            CalcH8K(ref iInt, afE, faB, faWeight, faK);
            return dofEnumerator.GetTransformedMatrix(new SymmetricMatrix2D(faK));
        }

        public IMatrix2D MassMatrix(IElement element)
        {
            double[,] faXYZ = GetCoordinates(element);
            double[,] faDS = new double[iInt3, 24];
            double[,] faS = new double[iInt3, 8];
            double[, ,] faB = new double[iInt3, 24, 6];
            double[] faDetJ = new double[iInt3];
            double[, ,] faJ = new double[iInt3, 3, 3];
            double[] faWeight = new double[iInt3];
            double[] faM = new double[300];
            double fDensity = Density;
            CalcH8GaussMatrices(ref iInt, faXYZ, faWeight, faS, faDS, faJ, faDetJ, faB);
            CalcH8MLumped(ref iInt, ref fDensity, faWeight, faM);
            return new SymmetricMatrix2D(faM);
        }

        public IMatrix2D DampingMatrix(IElement element)
        {
            var m = MassMatrix(element);
            var lc = m as ILinearlyCombinable;
            lc.LinearCombination(new double[] { RayleighAlpha, RayleighBeta }, new IMatrix2D[] { MassMatrix(element), StiffnessMatrix(element) });
            return m;
        }

        private static int GetSolidDOFs()
        {
            int totalDisplacementDOFs = 0;
            foreach (DOFType[] nodalDOFs in dofTypes)
                foreach (DOFType dofType in nodalDOFs)
                    if (dofType != DOFType.Pore) totalDisplacementDOFs++;
            return totalDisplacementDOFs;
        }

        private static int GetAllDOFs()
        {
            int totalDisplacementDOFs = 0;
            foreach (DOFType[] nodalDOFs in dofTypes)
                foreach (DOFType dofType in nodalDOFs)
                    totalDisplacementDOFs++;
            return totalDisplacementDOFs;
        }

        private static double[] ExtractSolidVector(double[] localVector)
        {
            int localPos = 0;
            int solidPos = 0;
            double[] solidVector = new double[GetSolidDOFs()];
            foreach (DOFType[] nodalDOFs in dofTypes)
                foreach (DOFType dofType in nodalDOFs)
                {
                    if (dofType != DOFType.Pore)
                    {
                        solidVector[solidPos] = localVector[localPos];
                        solidPos++;
                    }
                    localPos++;
                }

            return solidVector;
        }

        private static double[] ExtractFluidVector(double[] localVector)
        {
            int localPos = 0;
            int fluidPos = 0;
            double[] fluidVector = new double[GetAllDOFs() - GetSolidDOFs()];
            foreach (DOFType[] nodalDOFs in dofTypes)
                foreach (DOFType dofType in nodalDOFs)
                {
                    if (dofType == DOFType.Pore)
                    {
                        fluidVector[fluidPos] = localVector[localPos];
                        fluidPos++;
                    }
                    localPos++;
                }

            return fluidVector;
        }

        private static void ScatterFromSolidVector(double[] solidVector, double[] totalVector)
        {
            int localPos = 0;
            int solidPos = 0;
            foreach (DOFType[] nodalDOFs in dofTypes)
                foreach (DOFType dofType in nodalDOFs)
                {
                    if (dofType != DOFType.Pore)
                    {
                        totalVector[localPos] = solidVector[solidPos];
                        solidPos++;
                    }
                    localPos++;
                }
        }

        private static void ScatterFromFluidVector(double[] fluidVector, double[] totalVector)
        {
            int localPos = 0;
            int fluidPos = 0;
            foreach (DOFType[] nodalDOFs in dofTypes)
                foreach (DOFType dofType in nodalDOFs)
                {
                    if (dofType == DOFType.Pore)
                    {
                        totalVector[localPos] = fluidVector[fluidPos];
                        fluidPos++;
                    }
                    localPos++;
                }
        }

        public Tuple<double[], double[]> CalculateStresses(Element element, double[] localDisplacements, double[] localdDisplacements)
        {
            double[,] faXYZ = GetCoordinates(element);
            double[,] faDS = new double[iInt3, 24];
            double[,] faS = new double[iInt3, 8];
            double[, ,] faB = new double[iInt3, 24, 6];
            double[] faDetJ = new double[iInt3];
            double[, ,] faJ = new double[iInt3, 3, 3];
            double[] faWeight = new double[iInt3];
            double[,] faStrains = new double[iInt3, 6];
            double[,] fadStrains = new double[iInt3, 6];
            CalcH8GaussMatrices(ref iInt, faXYZ, faWeight, faS, faDS, faJ, faDetJ, faB);

            double[] soliddDisplacements = ExtractSolidVector(localdDisplacements);
            double[] solidDisplacements = ExtractSolidVector(localDisplacements);
            CalcH8Strains(ref iInt, faB, soliddDisplacements, fadStrains);
            CalcH8Strains(ref iInt, faB, solidDisplacements, faStrains);

            double[] dStrains = new double[6];
            double[] strains = new double[6];
            for (int i = 0; i < materialsAtGaussPoints.Length; i++)
            {
                for (int j = 0; j < 6; j++) dStrains[j] = fadStrains[i, j];
                for (int j = 0; j < 6; j++) strains[j] = faStrains[i, j];
                materialsAtGaussPoints[i].UpdateMaterial(new StressStrainVectorContinuum3D(dStrains));
            }

            return new Tuple<double[], double[]>(strains, materialsAtGaussPoints[materialsAtGaussPoints.Length - 1].Stresses.Data);
        }

        public double[] CalculateForcesForLogging(Element element, double[] localDisplacements)
        {
            return CalculateForces(element, localDisplacements, new double[localDisplacements.Length]);
        }

        public double[] CalculateForces(Element element, double[] localTotalDisplacements, double[] localDisplacements)
        {
            double[,] faStresses = new double[iInt3, 6];
            for (int i = 0; i < materialsAtGaussPoints.Length; i++)
                for (int j = 0; j < 6; j++) faStresses[i, j] = materialsAtGaussPoints[i].Stresses[j];

            double[,] faXYZ = GetCoordinates(element);
            double[,] faDS = new double[iInt3, 24];
            double[,] faS = new double[iInt3, 8];
            double[, ,] faB = new double[iInt3, 24, 6];
            double[] faDetJ = new double[iInt3];
            double[, ,] faJ = new double[iInt3, 3, 3];
            double[] faWeight = new double[iInt3];
            double[] solidForces = new double[24];
            CalcH8GaussMatrices(ref iInt, faXYZ, faWeight, faS, faDS, faJ, faDetJ, faB);
            CalcH8Forces(ref iInt, faB, faWeight, faStresses, solidForces);

            double[] faH = new double[36];
            double[] faPermeability = new double[] { Permeability, Permeability, Permeability, Permeability, 
                Permeability, Permeability, Permeability, Permeability };
            CalcH20u8pH(ref iInt, faPermeability, faS, faB, faWeight, faH);

            //double[,] faQ = new double[8, 24];
            //double fPoreA = PoreA;
            //double[] faXw = new double[] { Xw, Xw, Xw, Xw, Xw, Xw, Xw, Xw };
            //CalcH8u8pQMinus(ref iInt, ref fPoreA, faXw, faB, faS, faWeight, faQ);
            //Matrix2D<double> q = (new Matrix2D<double>(faQ)).Transpose();

            // Changed! Check for errors...
            //Vector<double> fluidDisplacements = new Vector<double>(ExtractFluidVector(localDisplacements));
            Vector fluidDisplacements = new Vector(ExtractFluidVector(localTotalDisplacements));
            
            //double[] solidAndFluidForces = solidForces;
            ////double[] solidAndFluidForces = q * fluidDisplacements + (new Vector<double>(solidForces));
            // H correction
            fluidDisplacements.Scale(-1);
            double[] fluidDrags = (new SymmetricMatrix2D(faH)) * fluidDisplacements;

            double[] totalForces = new double[GetAllDOFs()];
            ScatterFromFluidVector(fluidDrags, totalForces);
            ScatterFromSolidVector(solidForces, totalForces);
            return totalForces;
        }

        public double[] CalculateAccelerationForces(Element element, IList<MassAccelerationLoad> loads)
        {
            Vector accelerations = new Vector(24);
            int index = 0;
            foreach (MassAccelerationLoad load in loads)
                switch (load.DOF)
                {
                    case DOFType.X:
                        for (int i = 0; i < 8; i++) accelerations[i * 3] += load.Amount;
                        break;
                    case DOFType.Y:
                        for (int i = 0; i < 8; i++) accelerations[i * 3 + 1] += load.Amount;
                        break;
                    case DOFType.Z:
                        for (int i = 0; i < 8; i++) accelerations[i * 3 + 2] += load.Amount;
                        break;
                    default:
                        throw new InvalidOperationException("Cannot handle global acceleration for water pore when NOT translational.");
                }
            double[] solidForces = new double[24];
            MassMatrix(element).Multiply(accelerations, solidForces);

            double[,] faXYZ = GetCoordinates(element);
            double[,] faDS = new double[iInt3, 24];
            double[,] faS = new double[iInt3, 8];
            double[, ,] faB = new double[iInt3, 24, 6];
            double[] faDetJ = new double[iInt3];
            double[, ,] faJ = new double[iInt3, 3, 3];
            double[] faWeight = new double[iInt3];
            CalcH8GaussMatrices(ref iInt, faXYZ, faWeight, faS, faDS, faJ, faDetJ, faB);

            double[] waterAcc = new double[3];
            foreach (MassAccelerationLoad load in loads)
                switch (load.DOF)
                {
                    case DOFType.X:
                        waterAcc[0] = load.Amount;
                        break;
                    case DOFType.Y:
                        waterAcc[1] = load.Amount;
                        break;
                    case DOFType.Z:
                        waterAcc[2] = load.Amount;
                        break;
                    default:
                        throw new InvalidOperationException("Cannot handle global acceleration for water pore when NOT translational.");
                }

            bool[] impermeableDOFs = new bool[8];
            index = 0;
            foreach (Node node in element.NodesDictionary.Values)
                foreach (DOFType dofType in element.Subdomain.NodalDOFsDictionary[node.ID].Keys)
                {
                    if (dofType != DOFType.Pore) continue;
                    if (element.Subdomain.NodalDOFsDictionary[node.ID][dofType] < 0) impermeableDOFs[index] = true;
                    index++;
                }

            double fD = FluidDensity;
            double[] faPermeability = new double[] { Permeability, Permeability, Permeability, Permeability, 
                Permeability, Permeability, Permeability, Permeability };
            double[] faXw = new double[] { Xw, Xw, Xw, Xw, Xw, Xw, Xw, Xw };
            double[] faSw = new double[] { Saturation, Saturation, Saturation, Saturation, Saturation, 
                Saturation, Saturation, Saturation };
            double[] fluidDrags = new double[8];
            CalcH20u8pForcesWaterAcc(ref iInt, impermeableDOFs, ref fD, faPermeability,
                faXw, faSw, waterAcc, faS, faB, faWeight, fluidDrags);

            // H correction
            for (int i = 0; i < 8; i++)
                fluidDrags[i] = -fluidDrags[i];
            double[] totalForces = new double[GetAllDOFs()];
            ScatterFromFluidVector(fluidDrags, totalForces);
            ScatterFromSolidVector(solidForces, totalForces);
            return totalForces;
        }

        public double[] CalculateSolidForcesFromPorePressures(Element element, double[] porePressures)
        {
            double[] solidForces = new double[24];
            Matrix2D Q = ((Matrix2D)CouplingMatrix(element)).Transpose();
            Q.Multiply(new Vector(porePressures), solidForces);
            return solidForces;
        }

        //public IFiniteElementMaterial Material
        //{
        //    get { return materialsAtGaussPoints; }
        //}

        public bool MaterialModified
        {
            get 
            {
                foreach (IContinuumMaterial3D material in materialsAtGaussPoints)
                    if (material.Modified) return true;
                return false;
            }
        }

        public void ResetMaterialModified()
        {
            foreach (IContinuumMaterial3D material in materialsAtGaussPoints) material.ResetModified();
        }

        public void ClearMaterialState()
        {
            foreach (IContinuumMaterial3D m in materialsAtGaussPoints) m.ClearState();
        }

        public void SaveMaterialState()
        {
            foreach (IContinuumMaterial3D m in materialsAtGaussPoints) m.SaveState();
        }

        public void ClearMaterialStresses()
        {
            foreach (IContinuumMaterial3D m in materialsAtGaussPoints) m.ClearStresses();
        }

        #endregion

        #region IPorousFiniteElement Members

        public IMatrix2D PermeabilityMatrix(IElement element)
        {
            double[,] faXYZ = GetCoordinates(element);
            double[,] faDS = new double[iInt3, 24];
            double[,] faS = new double[iInt3, 8];
            double[, ,] faB = new double[iInt3, 24, 6];
            double[] faDetJ = new double[iInt3];
            double[, ,] faJ = new double[iInt3, 3, 3];
            double[] faWeight = new double[iInt3];
            double[] faPermeability = new double[] { Permeability, Permeability, Permeability, Permeability, 
                Permeability, Permeability, Permeability, Permeability };
            double[] faH = new double[36];
            CalcH8GaussMatrices(ref iInt, faXYZ, faWeight, faS, faDS, faJ, faDetJ, faB);
            CalcH20u8pH(ref iInt, faPermeability, faS, faB, faWeight, faH);
            return new SymmetricMatrix2D(faH);
        }

        // Rows are fluid DOFs and columns are solid DOFs
        public IMatrix2D CouplingMatrix(IElement element)
        {
            double[,] faXYZ = GetCoordinates(element);
            double[,] faDS = new double[iInt3, 24];
            double[,] faS = new double[iInt3, 8];
            double[, ,] faB = new double[iInt3, 24, 6];
            double[] faDetJ = new double[iInt3];
            double[, ,] faJ = new double[iInt3, 3, 3];
            double[] faWeight = new double[iInt3];
            //double fDensity = Density;
            double fPoreA = PoreA;
            double[] faXw = new double[] { Xw, Xw, Xw, Xw, Xw, Xw, Xw, Xw };
            double[,] faQ = new double[8, 24];
            CalcH8GaussMatrices(ref iInt, faXYZ, faWeight, faS, faDS, faJ, faDetJ, faB);
            CalcH8u8pQMinus(ref iInt, ref fPoreA, faXw, faB, faS, faWeight, faQ);
            return new Matrix2D(faQ);
        }

        public IMatrix2D SaturationMatrix(IElement element)
        {
            double[,] faXYZ = GetCoordinates(element);
            double[,] faDS = new double[iInt3, 24];
            double[,] faS = new double[iInt3, 8];
            double[, ,] faB = new double[iInt3, 24, 6];
            double[] faDetJ = new double[iInt3];
            double[, ,] faJ = new double[iInt3, 3, 3];
            double[] faWeight = new double[iInt3];
            double[] faXwDivQ = new double[] { Xw * QInv, Xw * QInv, Xw * QInv, Xw * QInv, Xw * QInv, 
                Xw * QInv, Xw * QInv, Xw * QInv};
            double[] faSaturation = new double[36];
            CalcH8GaussMatrices(ref iInt, faXYZ, faWeight, faS, faDS, faJ, faDetJ, faB);
            CalcH20u8pS(ref iInt, faXwDivQ, faS, faWeight, faSaturation);
            return new SymmetricMatrix2D(faSaturation);
        }

        #endregion
    }
}
