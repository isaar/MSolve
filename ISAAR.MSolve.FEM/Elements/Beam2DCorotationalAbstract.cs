using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Elements.SupportiveClasses;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.Materials.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;

namespace ISAAR.MSolve.FEM.Elements
{
    public abstract class Beam2DCorotationalAbstract : IFiniteElement
    {
        protected static readonly int NATURAL_DEFORMATION_COUNT = 3;
        protected static readonly int FREEDOM_DEGREE_COUNT = 6;
        protected static readonly int AXIS_COUNT = 1;
        protected static readonly int NODE_COUNT = 2;

        protected readonly static DOFType[] nodalDOFTypes = new DOFType[] { DOFType.X, DOFType.Y, DOFType.RotZ };
        protected readonly static DOFType[][] dofTypes = new DOFType[][] { nodalDOFTypes, nodalDOFTypes };
        private static readonly DOFType[][] dofs = new DOFType[][] { nodalDOFTypes, nodalDOFTypes };

        protected IElementDOFEnumerator dofEnumerator = new GenericDOFEnumerator();
        protected readonly IFiniteElementMaterial material;
        protected readonly IList<Node> nodes;
        protected readonly double density;
        protected BeamSection2D beamSection;
        protected readonly double initialLength;
        protected double currentLength;
        protected Matrix2D currentRotationMatrix;
        protected Vector naturalDeformations;
        protected Vector beamAxisX;
        protected Vector beamAxisY;
        protected Vector beamAxisZ;

        public double RayleighAlpha { get; set; }
        public double RayleighBeta { get; set; }

        protected Beam2DCorotationalAbstract(IList<Node> nodes, IFiniteElementMaterial material, double density, BeamSection2D beamSection)
        {
            this.nodes = nodes;
            this.material = material;
            this.density = density;
            this.beamSection = beamSection;
            this.initialLength = Math.Sqrt(Math.Pow(nodes[0].X - nodes[1].X, 2) + Math.Pow(nodes[0].Y - nodes[1].Y, 2));
            this.currentLength = this.initialLength;
            this.currentRotationMatrix = new Matrix2D(AXIS_COUNT, AXIS_COUNT);
            this.naturalDeformations = new Vector(NATURAL_DEFORMATION_COUNT);
            this.beamAxisX = new Vector(AXIS_COUNT);
            this.beamAxisY = new Vector(AXIS_COUNT);
        }

        public int ID { get { return 100; } }
        public ElementDimensions ElementDimensions { get { return ElementDimensions.ThreeD; } }
        public bool MaterialModified { get { return material.Modified; } }
        public IElementDOFEnumerator DOFEnumerator
        {
            get { return dofEnumerator; }
            set { dofEnumerator = value; }
        }

        public abstract void SaveGeometryState();
        public abstract void UpdateState(double[] incrementalNodeDisplacements);

        private Matrix2D CalculateBlockRotationMatrix()
        {
            Matrix2D blockRotationMatrix = new Matrix2D(FREEDOM_DEGREE_COUNT, FREEDOM_DEGREE_COUNT);
            int totalBlocks = 2;
            int blockSize = 3;
            double R11 = this.currentRotationMatrix[0, 0];
            double R12 = this.currentRotationMatrix[0, 1];
            double R13 = this.currentRotationMatrix[0, 2];
            double R21 = this.currentRotationMatrix[1, 0];
            double R22 = this.currentRotationMatrix[1, 1];
            double R23 = this.currentRotationMatrix[1, 2];
            double R31 = this.currentRotationMatrix[2, 0];
            double R32 = this.currentRotationMatrix[2, 1];
            double R33 = this.currentRotationMatrix[2, 2];

            for (int block = 0; block < totalBlocks; block++)
            {
                int cellDistanceCount = block * blockSize;

                blockRotationMatrix[cellDistanceCount, cellDistanceCount] = R11;
                blockRotationMatrix[cellDistanceCount, cellDistanceCount + 1] = R12;
                blockRotationMatrix[cellDistanceCount, cellDistanceCount + 2] = R13;

                blockRotationMatrix[cellDistanceCount + 1, cellDistanceCount] = R21;
                blockRotationMatrix[cellDistanceCount + 1, cellDistanceCount + 1] = R22;
                blockRotationMatrix[cellDistanceCount + 1, cellDistanceCount + 2] = R23;

                blockRotationMatrix[cellDistanceCount + 2, cellDistanceCount] = R31;
                blockRotationMatrix[cellDistanceCount + 2, cellDistanceCount + 1] = R32;
                blockRotationMatrix[cellDistanceCount + 2, cellDistanceCount + 2] = R33;
            }

            return blockRotationMatrix;
        }

        /**
	     * Calculates the constitutive stiffness of the element.
	     *
	     * @return The constitutive stiffness
	     */
        private SymmetricMatrix2D CalculateConstitutiveStiffness()
        {
            var constitutiveStiffness = new SymmetricMatrix2D(FREEDOM_DEGREE_COUNT);
            double E = this.material.YoungModulus;
            double G = E / (2d * (1d + this.material.PoissonRatio));
            double I = this.beamSection.Inertia;
            double A = this.beamSection.Area;
            double L = this.currentLength;
            double LSqruared = L * L;
            double LCubed = L * L * L;
            double phi = (12.0 * E * I) / (LSqruared * G * A);
            double psi = 1.0 / (1.0 + phi);
            double EAOverL = (E * A) / L;
            double EIOverL = (E * I) / L;

            constitutiveStiffness[0, 0] = EAOverL;
            constitutiveStiffness[0, 3] = -EAOverL;

            constitutiveStiffness[1, 1] = 12.0 * psi * E * I / LCubed;
            constitutiveStiffness[1, 2] = 6.0 * psi * E * I / LSqruared;
            constitutiveStiffness[1, 4] = -12.0 * psi * E * I / LCubed;
            constitutiveStiffness[1, 5] = 6.0 * psi * E * I / LSqruared;

            constitutiveStiffness[2, 2] = (3.0 * psi + 1.0) * EIOverL;
            constitutiveStiffness[2, 4] = -6.0 * psi * E * I / L;
            constitutiveStiffness[2, 5] = (3.0 * psi + 1.0) * EIOverL;

            constitutiveStiffness[3, 3] = EAOverL;

            constitutiveStiffness[4, 4] = 12.0 * psi * E * I / LCubed;
            constitutiveStiffness[4, 5] = -6.0 * psi * E * I / LSqruared;

            constitutiveStiffness[5, 5] = (3.0 * psi + 1.0) * EIOverL;

            return constitutiveStiffness;
        }

        /**
	     * Calculates the forces in the global coordinate system.
	     *
	     * @return The forces in the global coordinate system
	     */
        private Vector CalculateForcesInGlobalSystem()
        {
            var forcesNatural = this.CalculateForcesInNaturalSystem();
            var transformationMatrix = this.CalculateNaturalToGlobalTransormMatrix();
            double[] forcesGlobal = new double[FREEDOM_DEGREE_COUNT];
            transformationMatrix.Multiply(forcesNatural, forcesGlobal);
            return new Vector(forcesGlobal);
        }

        /**
	     * Calculates the forces in the local coordinate system.
	     *
	     * @return The forces in the local coordinate system
	     */
        private Vector CalculateForcesInLocalSystem()
        {
            var naturalToLocal = this.CalculateNaturalToLocalTranformMatrix();
            var naturalForces = this.CalculateForcesInNaturalSystem();
            double[] forcesLocal = new double[FREEDOM_DEGREE_COUNT];
            naturalToLocal.Multiply(naturalForces, forcesLocal);
            return new Vector(forcesLocal);
        }

        /**
	     * Calculates forces in the natural coordinate system.
	     *
	     * @return The forces in the natural coordinate system
	     */
        private Vector CalculateForcesInNaturalSystem()
        {
            var forcesNatural = new Vector(NATURAL_DEFORMATION_COUNT);
            double E = this.material.YoungModulus;
            double G = E / (2d * (1d + this.material.PoissonRatio));
            double I = this.beamSection.Inertia;
            double A = this.beamSection.Area;
            double l = this.currentLength;
            double phi = (12.0 * E * I) / (l * l * G * A);
            double psi = 1.0 / (1.0 + phi);

            forcesNatural[NaturalDeformationMode2D.EXTENSION] =
                (E * A * this.naturalDeformations[NaturalDeformationMode2D.EXTENSION]) / l;

            forcesNatural[NaturalDeformationMode2D.SYMMETRIC_BENDING] =
                (E * I * this.naturalDeformations[NaturalDeformationMode2D.SYMMETRIC_BENDING]) / l;
            
            forcesNatural[NaturalDeformationMode2D.ANTISYMMETRIC_BENDING] =
                (3.0 * psi * E * I * this.naturalDeformations[NaturalDeformationMode2D.ANTISYMMETRIC_BENDING]) / l;

            return forcesNatural;
        }

        /**
	     * Calculates the geometric stiffness of the element.
	     *
	     * @return The geometric stiffness
	     */
        private SymmetricMatrix2D CalculateGeometricStiffness()
        {
            var forcesNatural = new Vector(NATURAL_DEFORMATION_COUNT);
            double E = this.material.YoungModulus;
            double G = E / (2d * (1d + this.material.PoissonRatio));
            double I = this.beamSection.Inertia;
            double A = this.beamSection.Area;
            double L = this.currentLength;
            var geometricStiffness = new SymmetricMatrix2D(FREEDOM_DEGREE_COUNT);
            var forcesInNaturalSystem = this.CalculateForcesInNaturalSystem();
            double axialForce = forcesInNaturalSystem[NaturalDeformationMode2D.EXTENSION];
            double shearForce = -2.0*forcesInNaturalSystem[2] / L;

            geometricStiffness[0, 1] = -shearForce / L;
            geometricStiffness[0, 4] = shearForce / L;

            geometricStiffness[1, 1] = 6.0 * axialForce / (5.0 * L);
            geometricStiffness[1, 2] = axialForce / 10.0;
            geometricStiffness[1, 3] = shearForce / L;
            geometricStiffness[1, 4] = -6.0 * axialForce / (5.0 * L);
            geometricStiffness[1, 5] = axialForce / 10.0;

            geometricStiffness[2, 2] = 2.0 * axialForce * L / 15.0;
            geometricStiffness[2, 4] = -axialForce / 10.0;
            geometricStiffness[2, 5] = -axialForce * L / 30.0;

            geometricStiffness[3, 4] = -shearForce / L;
            geometricStiffness[4, 4] = 6.0 * axialForce / (5.0 * L);
            geometricStiffness[4, 5] = -axialForce / 10.0;
            geometricStiffness[5, 5] = 2.0 * axialForce * L / 15.0;
            
            return geometricStiffness;
        }

        /**
	     * Calculates the stiffness matrix in the local coordinate system.
	     *
	     * @return The stiffness matrix in the local coordinate system.
	     */
        private SymmetricMatrix2D CalculateLocalStiffnessMatrix()
        {
            var constitutivePart = this.CalculateConstitutiveStiffness();
            var geometricPart = this.CalculateGeometricStiffness();
            constitutivePart.LinearCombination(new[] { 1d, 1d }, new List<SymmetricMatrix2D> { constitutivePart, geometricPart });
            return constitutivePart;
        }

        /**
	     * Calculates the transformation matrix from natural to local coordinate system.
	     *
	     * @return The natural to local transformation matrix
	     */
        private Matrix2D CalculateNaturalToGlobalTransormMatrix()
        {
            var transformMatrix = new Matrix2D(FREEDOM_DEGREE_COUNT, NATURAL_DEFORMATION_COUNT);
            double L = this.currentLength;

            transformMatrix[0, 0] = -1.0;
            transformMatrix[1, 2] = +2.0 / L;
            transformMatrix[2, 1] = -1.0;
            transformMatrix[2, 2] = +1.0;
            transformMatrix[3, 0] = +1.0;
            transformMatrix[4, 2] = -2.0 / L;
            transformMatrix[5, 1] = +1.0;
            transformMatrix[5, 2] = +1.0;
            
            return transformMatrix;
        }

        /**
	     * Calculates the transformation matrix from natural to local coordinate system.
	     *
	     * @return The natural to local transformation matrix
	     */
        private Matrix2D CalculateNaturalToLocalTranformMatrix()
        {
            var transformMatrix = new Matrix2D(FREEDOM_DEGREE_COUNT, NATURAL_DEFORMATION_COUNT);
            double L = this.currentLength;

            transformMatrix[0, 0] = -1.0;
            transformMatrix[1, 2] = +2.0 / L;
            transformMatrix[2, 1] = -1.0;
            transformMatrix[2, 2] = +1.0;
            transformMatrix[3, 0] = +1.0;
            transformMatrix[4, 2] = -2.0 / L;
            transformMatrix[5, 1] = +1.0;
            transformMatrix[5, 2] = +1.0;

            return transformMatrix;
        }

        public IList<IList<DOFType>> GetElementDOFTypes(IElement element)
        {
            return dofTypes;
        }

        public IMatrix2D StiffnessMatrix(IElement element)
        {
            var rotationMatrixBlock = this.CalculateBlockRotationMatrix();
            var localStiffnessMatrix = this.CalculateLocalStiffnessMatrix();
            var s = rotationMatrixBlock * localStiffnessMatrix.ToMatrix2D() * rotationMatrixBlock.Transpose();
            return new SymmetricMatrix2D(s);
        }

        public IMatrix2D MassMatrix(IElement element)
        {
            throw new NotImplementedException();
            //double area = beamSection.Area;
            //double inertiaY = beamSection.InertiaY;
            //double inertiaZ = beamSection.InertiaZ;
            //double x2 = Math.Pow(element.INodes[1].X - element.INodes[0].X, 2);
            //double y2 = Math.Pow(element.INodes[1].Y - element.INodes[0].Y, 2);
            //double z2 = Math.Pow(element.INodes[1].Z - element.INodes[0].Z, 2);
            //double L = Math.Sqrt(x2 + y2 + z2);
            //double fullMass = density * area * L;

            //var massMatrix = new Matrix2D(FREEDOM_DEGREE_COUNT, FREEDOM_DEGREE_COUNT);
            //massMatrix[0, 0] = (1.0 / 3.0) * fullMass;
            //massMatrix[0, 6] = (1.0 / 6.0) * fullMass;

            //massMatrix[1, 1] = (13.0 / 35.0) * fullMass;
            //massMatrix[1, 5] = (11.0 * L / 210.0) * fullMass;
            //massMatrix[1, 7] = (9.0 / 70.0) * fullMass;
            //massMatrix[1, 11] = -(13.0 * L / 420.0) * fullMass;

            //massMatrix[2, 2] = (13.0 / 35.0) * fullMass;
            //massMatrix[2, 4] = -(11.0 * L / 210.0) * fullMass;
            //massMatrix[2, 8] = (9.0 / 70.0) * fullMass;
            //massMatrix[2, 10] = -(13.0 * L / 420.0) * fullMass;

            //massMatrix[3, 3] = ((inertiaY + inertiaZ) / (3.0 * area)) * fullMass;
            //massMatrix[3, 9] = ((inertiaY + inertiaZ) / (6.0 * area)) * fullMass;

            //massMatrix[4, 4] = ((L * L) / 105.0) * fullMass;
            //massMatrix[4, 8] = -(13.0 * L / 420.0) * fullMass;
            //massMatrix[4, 10] = -((L * L) / 105.0) * fullMass;

            //massMatrix[5, 5] = ((L * L) / 105.0) * fullMass;
            //massMatrix[5, 7] = (13.0 * L / 420.0) * fullMass;
            //massMatrix[5, 11] = -((L * L) / 105.0) * fullMass;

            //massMatrix[6, 6] = (1.0 / 3.0) * fullMass;

            //massMatrix[7, 7] = (13.0 / 35.0) * fullMass;
            //massMatrix[7, 11] = -(11.0 * L / 210.0) * fullMass;

            //massMatrix[8, 8] = (13.0 / 35.0) * fullMass;
            //massMatrix[8, 10] = (11.0 * L / 210.0) * fullMass;

            //massMatrix[9, 9] = ((inertiaY + inertiaZ) / (3.0 * area)) * fullMass;

            //massMatrix[10, 10] = ((L * L) / 105.0) * fullMass;

            //massMatrix[11, 11] = ((L * L) / 105.0) * fullMass;

            //massMatrix[6, 0] = (1.0 / 6.0) * fullMass;
            //massMatrix[5, 1] = (11.0 * L / 210.0) * fullMass;
            //massMatrix[7, 1] = (9.0 / 70.0) * fullMass;
            //massMatrix[11, 1] = -(13.0 * L / 420.0) * fullMass;
            //massMatrix[4, 2] = -(11.0 * L / 210.0) * fullMass;
            //massMatrix[8, 2] = (9.0 / 70.0) * fullMass;
            //massMatrix[10, 2] = -(13.0 * L / 420.0) * fullMass;
            //massMatrix[9, 3] = ((inertiaY + inertiaZ) / (6.0 * area)) * fullMass;
            //massMatrix[8, 4] = -(13.0 * L / 420.0) * fullMass;
            //massMatrix[10, 4] = -((L * L) / 105.0) * fullMass;
            //massMatrix[7, 5] = (13.0 * L / 420.0) * fullMass;
            //massMatrix[11, 5] = -((L * L) / 105.0) * fullMass;
            //massMatrix[11, 7] = -(11.0 * L / 210.0) * fullMass;
            //massMatrix[10, 8] = (11.0 * L / 210.0) * fullMass;

            //return massMatrix;
        }

        public IMatrix2D DampingMatrix(IElement element)
        {
            throw new NotImplementedException();
            //var m = MassMatrix(element);
            //var lc = m as ILinearlyCombinable;
            //lc.LinearCombination(new double[] { RayleighAlpha, RayleighBeta }, new IMatrix2D[] { MassMatrix(element), StiffnessMatrix(element) });
            //return m;
        }

        public void ResetMaterialModified()
        {
            this.material.ResetModified();
        }

        public Tuple<double[], double[]> CalculateStresses(Element element, double[] localDisplacements, double[] localdDisplacements)
        {
            UpdateState(localdDisplacements);
            //TODO: Should calculate strains and update material as well
            //material.UpdateMaterial(strains);
            //TODO: Should calculate stresses as well
            return new Tuple<double[], double[]>(new double[FREEDOM_DEGREE_COUNT], new double[FREEDOM_DEGREE_COUNT]);
        }

        public double[] CalculateForces(Element element, double[] localDisplacements, double[] localdDisplacements)
        {
            var internalForces = this.CalculateForcesInGlobalSystem();
            return internalForces.Data;
        }

        public double[] CalculateForcesForLogging(Element element, double[] localDisplacements)
        {
            throw new NotImplementedException();
        }

        public double[] CalculateAccelerationForces(Element element, IList<MassAccelerationLoad> loads)
        {
            Vector accelerations = new Vector(6);
            IMatrix2D massMatrix = MassMatrix(element);

            int index = 0;
            foreach (MassAccelerationLoad load in loads)
                foreach (DOFType[] nodalDOFTypes in dofs)
                    foreach (DOFType dofType in nodalDOFTypes)
                    {
                        if (dofType == load.DOF) accelerations[index] += load.Amount;
                        index++;
                    }

            double[] forces = new double[6];
            massMatrix.Multiply(accelerations, forces);
            return forces;
        }

        public void SaveMaterialState()
        {
            SaveGeometryState();
            material.SaveState();
        }

        public void ClearMaterialState()
        {
            material.ClearState();
        }

        public void ClearMaterialStresses()
        {
            material.ClearStresses();
        }
    }
}
