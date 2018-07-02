using ISAAR.MSolve.FEM.Elements.SupportiveClasses;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.Materials.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.Matrices;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.FEM.Elements
{
    public class Beam3DCorotationalQuaternion : Beam3DCorotationalAbstract
    {
        private readonly Vector lastDisplacements;
    	private readonly Vector currentDisplacements;
        private readonly Vector displacementsOfCurrentIncrement;
        private readonly Vector initialAxisY;
    	private readonly Vector initialAxisZ;
    	private readonly Vector initialAxisX;
    	private Quaternion quaternionCurrentNodeA;
        private Quaternion quaternionCurrentNodeB;
        private readonly Vector currentBeamAxis;
        private Quaternion quaternionLastNodeA;
        private Quaternion quaternionLastNodeB;

        /**
         * Creates a new instance of {@link Beam3DCorotationalIncremental} class.
         *
         * @param nodes
         *            The element nodes
         * @param material
         *            The element material
         * @param density
         *            The element density
         * @param beamSection
         *            The beam section.
         */
        public Beam3DCorotationalQuaternion(IList<Node> nodes, IIsotropicContinuumMaterial3D material, double density, BeamSection3D beamSection)
            : base(nodes, material, density, beamSection)
        {
            this.displacementsOfCurrentIncrement = new Vector(FREEDOM_DEGREE_COUNT);
            this.lastDisplacements = new Vector(FREEDOM_DEGREE_COUNT);
            this.currentDisplacements = new Vector(FREEDOM_DEGREE_COUNT);
            this.quaternionLastNodeA = Quaternion.OfZeroAngle();
            this.quaternionLastNodeB = Quaternion.OfZeroAngle();
            this.quaternionCurrentNodeA = Quaternion.OfZeroAngle();
            this.quaternionCurrentNodeB = Quaternion.OfZeroAngle();
            this.initialAxisX = new Vector(AXIS_COUNT);
            this.initialAxisY = new Vector(AXIS_COUNT);
            this.initialAxisZ = new Vector(AXIS_COUNT);
            this.currentBeamAxis = new Vector(AXIS_COUNT);
            this.InitializeElementAxes();
        }

        public override void SaveGeometryState()
        {
            displacementsOfCurrentIncrement.Scale(0d);
            currentDisplacements.CopyTo(lastDisplacements.Data, 0);

            var qA = new double[quaternionCurrentNodeA.VectorPart.Length];
            quaternionCurrentNodeA.VectorPart.CopyTo(qA, 0);
            var qB = new double[quaternionCurrentNodeB.VectorPart.Length];
            quaternionCurrentNodeB.VectorPart.CopyTo(qB, 0);
            this.quaternionLastNodeA = new Quaternion(quaternionCurrentNodeA.ScalarPart, new Vector(qA));
            this.quaternionLastNodeB = new Quaternion(quaternionCurrentNodeB.ScalarPart, new Vector(qB));
        }

        public override void UpdateState(double[] incrementalNodeDisplacements)
        {
            //displacementsOfCurrentIncrement.Add(new Vector(incrementalNodeDisplacements));
            //lastDisplacements.CopyTo(currentDisplacements.Data, 0);
            //currentDisplacements.Add(displacementsOfCurrentIncrement);

            //this.currentDisplacements.addIntoThis(this.lastDisplacements, incrementalNodeDisplacements);
            lastDisplacements.CopyTo(currentDisplacements.Data, 0);
            currentDisplacements.Add(new Vector(incrementalNodeDisplacements));
            new Vector(incrementalNodeDisplacements).CopyTo(displacementsOfCurrentIncrement.Data, 0);

            var incrementalRotationsA = new Vector(AXIS_COUNT);
            var incrementalRotationsB = new Vector(AXIS_COUNT);

            incrementalRotationsA[0] = displacementsOfCurrentIncrement[3];
            incrementalRotationsA[1] = displacementsOfCurrentIncrement[4];
            incrementalRotationsA[2] = displacementsOfCurrentIncrement[5];
            incrementalRotationsB[0] = displacementsOfCurrentIncrement[9];
            incrementalRotationsB[1] = displacementsOfCurrentIncrement[10];
            incrementalRotationsB[2] = displacementsOfCurrentIncrement[11];

            var qA = new double[quaternionLastNodeA.VectorPart.Length];
            quaternionLastNodeA.VectorPart.CopyTo(qA, 0);
            var qB = new double[quaternionLastNodeB.VectorPart.Length];
            quaternionLastNodeB.VectorPart.CopyTo(qB, 0);
            this.quaternionCurrentNodeA = new Quaternion(quaternionLastNodeA.ScalarPart, new Vector(qA));
            this.quaternionCurrentNodeB = new Quaternion(quaternionLastNodeB.ScalarPart, new Vector(qB));
            this.quaternionCurrentNodeA.ApplyIncrementalRotation(incrementalRotationsA);
            this.quaternionCurrentNodeB.ApplyIncrementalRotation(incrementalRotationsB);

            double scalarA = this.quaternionCurrentNodeA.ScalarPart;
            double scalarB = this.quaternionCurrentNodeB.ScalarPart;
            var vectorPartA = this.quaternionCurrentNodeA.VectorPart;
            var vectorPartB = this.quaternionCurrentNodeB.VectorPart;
            var sumOfVectorParts = new Vector(vectorPartA.Length);
            vectorPartA.CopyTo(sumOfVectorParts.Data, 0);
            sumOfVectorParts.Add(vectorPartB);
            double sumOfVectorPartsNorm = sumOfVectorParts.Norm;
            double scalarPartDifference = 0.5 * Math.Sqrt(((scalarA + scalarB) * (scalarA + scalarB)) + (sumOfVectorPartsNorm * sumOfVectorPartsNorm));
            double meanRotationScalarPart = (0.5 * (scalarA + scalarB)) / scalarPartDifference;
            var meanRotationVectorPart = new Vector(vectorPartA.Length);
            vectorPartA.CopyTo(meanRotationVectorPart.Data, 0);
            meanRotationVectorPart.Add(vectorPartB);
            meanRotationVectorPart.Scale(1d / (2d * scalarPartDifference));
            var vectorPartDifference = new Vector(vectorPartB.Length);
            vectorPartB.CopyTo(vectorPartDifference.Data, 0);
            vectorPartDifference.Scale(scalarA);
            vectorPartDifference.Add(-scalarB * vectorPartA);
            //vectorPartDifference.doPointwise(vectorPartA, DoubleBinaryOps.alphaPlusScaledBeta(-scalarB));
            vectorPartDifference.Add(vectorPartA ^ vectorPartB);
            vectorPartDifference.Scale(1d / (2d * scalarPartDifference));
            Quaternion meanRotationQuaternion = Quaternion.CreateFromIndividualParts(meanRotationScalarPart, meanRotationVectorPart);
            this.currentRotationMatrix = meanRotationQuaternion.GetRotationMatrix();
            this.CalculateUpdatedBeamAxis();
            this.UpdateRotationMatrix();
            this.UpdateNaturalDeformations(vectorPartDifference);

            //SaveGeometryState();
        }

        private void CalculateUpdatedBeamAxis()
        {
            currentRotationMatrix.Multiply(initialAxisX, beamAxisX.Data);
            currentRotationMatrix.Multiply(initialAxisY, beamAxisY.Data);
            currentRotationMatrix.Multiply(initialAxisZ, beamAxisZ.Data);

            double dX = ((nodes[1].X - nodes[0].X) + this.currentDisplacements[6]) - this.currentDisplacements[0];
            double dY = ((nodes[1].Y - nodes[0].Y) + this.currentDisplacements[7]) - this.currentDisplacements[1];
            double dZ = ((nodes[1].Z - nodes[0].Z) + this.currentDisplacements[8]) - this.currentDisplacements[2];
            this.currentLength = Math.Sqrt((dX * dX) + (dY * dY) + (dZ * dZ));
            var delta = new Vector(new[] { dX, dY, dZ });
            this.currentBeamAxis[0] = dX / currentLength + beamAxisX[0];
            this.currentBeamAxis[1] = dY / currentLength + beamAxisX[1];
            this.currentBeamAxis[2] = dZ / currentLength + beamAxisX[2];
            this.currentBeamAxis.Scale(1d / this.currentBeamAxis.Norm);

            //this.currentBeamAxis.divideAllIntoThis(delta, this.currentLength);
            //this.currentBeamAxis.add(this.beamAxisX);
            //this.currentBeamAxis.divideAll(this.currentBeamAxis.normEntrywise(EntrywiseNorms.L2));
        }

        private void InitializeElementAxes()
        {
            var globalVectorX = new Vector(AXIS_COUNT);
            var globalVectorY = new Vector(AXIS_COUNT);
            var globalVectorZ = new Vector(AXIS_COUNT);
            globalVectorX[0] = 1d;
            globalVectorY[1] = 1d;
            globalVectorZ[2] = 1d;
            double deltaX = nodes[1].X - nodes[0].X;
            double deltaY = nodes[1].Y - nodes[0].Y;
            double deltaZ = nodes[1].Z - nodes[0].Z;
            var elementAxis = new Vector(3);
            elementAxis[0] = deltaX;
            elementAxis[1] = deltaY;
            elementAxis[2] = deltaZ;
            elementAxis.Scale(1d / elementAxis.Norm);

            currentRotationMatrix = RotationMatrix.CalculateRotationMatrix(globalVectorX, elementAxis);
            currentRotationMatrix.Multiply(globalVectorX, initialAxisX.Data);
            currentRotationMatrix.Multiply(globalVectorY, initialAxisY.Data);
            currentRotationMatrix.Multiply(globalVectorZ, initialAxisZ.Data);
            this.initialAxisX.CopyTo(beamAxisX.Data, 0);
            this.initialAxisY.CopyTo(beamAxisY.Data, 0);
            this.initialAxisZ.CopyTo(beamAxisZ.Data, 0);
        }

        private void UpdateNaturalDeformations(Vector vectorPartDifference)
        {
            double extension = this.currentLength - this.initialLength;
            Matrix2D currentRotationMatrixTranspose = currentRotationMatrix.Transpose();   
            var symmetricRotation = new Vector(AXIS_COUNT);

            //currentRotationMatrix.Transpose().Multiply(vectorPartDifference, symmetricRotation.Data);
            //currentRotationMatrix.Multiply(vectorPartDifference, symmetricRotation.Data);
            currentRotationMatrixTranspose.Multiply(vectorPartDifference, symmetricRotation.Data);

            var antisymmetricRotation = new Vector(AXIS_COUNT);

            //currentRotationMatrix.Multiply(this.beamAxisX ^ this.currentBeamAxis, antisymmetricRotation.Data);
            //currentRotationMatrix.Transpose().Multiply(this.beamAxisX ^ this.currentBeamAxis, antisymmetricRotation.Data);
            currentRotationMatrixTranspose.Multiply(this.beamAxisX ^ this.currentBeamAxis, antisymmetricRotation.Data);
            //currentRotationMatrix.MultiplyTranspose2(this.beamAxisX ^ this.currentBeamAxis, antisymmetricRotation.Data);
            //final VectorView symmetricRotation = vectorPartDifference.leftMultiplyWithMatrix(this.currentRotationMatrix);
            //final VectorView antisymmetricRotation =
            //        GeometricUtils.crossProduct(this.beamAxisX, this.currentBeamAxis)
            //                        .leftMultiplyWithMatrix(this.currentRotationMatrix);
            this.naturalDeformations[0] = symmetricRotation[0] * 4.0;
            this.naturalDeformations[1] = symmetricRotation[1] * 4.0;
            this.naturalDeformations[2] = symmetricRotation[2] * 4.0;
            this.naturalDeformations[3] = extension;
            this.naturalDeformations[4] = antisymmetricRotation[1] * 4.0;
            this.naturalDeformations[5] = antisymmetricRotation[2] * 4.0;
        }

        private void UpdateRotationMatrix()
        {
            var rotationMatrixLocal = new Matrix2D(AXIS_COUNT, AXIS_COUNT);
            for (int i = 0; i < AXIS_COUNT; i++)
                rotationMatrixLocal[i, i] = 1d;
            rotationMatrixLocal.LinearCombinationGOAT(new[] { -2d }, new[] { Matrix2D.FromVector(currentBeamAxis.Data) * Matrix2D.FromVectorTranspose(currentBeamAxis.Data) });
            var minusBeamAxisX = new Vector(beamAxisX.Length);
            var tempBeamAxisY = new Vector(beamAxisY.Length);
            var tempBeamAxisZ = new Vector(beamAxisZ.Length);
            beamAxisX.CopyTo(minusBeamAxisX.Data, 0);
            minusBeamAxisX.Scale(-1d);
            beamAxisY.CopyTo(tempBeamAxisY.Data, 0);
            beamAxisZ.CopyTo(tempBeamAxisZ.Data, 0);
            rotationMatrixLocal.Multiply(minusBeamAxisX, beamAxisX.Data);
            rotationMatrixLocal.Multiply(tempBeamAxisY, beamAxisY.Data);
            rotationMatrixLocal.Multiply(tempBeamAxisZ, beamAxisZ.Data);
            this.currentRotationMatrix = RotationMatrix.CalculateFromOrthonormalVectors(this.beamAxisX, this.beamAxisY, this.beamAxisZ);
        }

    }
}
