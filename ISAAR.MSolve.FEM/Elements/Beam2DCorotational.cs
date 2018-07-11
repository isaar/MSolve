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
    public class Beam2DCorotational : Beam2DCorotationalAbstract
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
        public Beam2DCorotational(IList<Node> nodes, IFiniteElementMaterial material, double density, BeamSection2D beamSection)
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
            lastDisplacements.CopyTo(currentDisplacements.Data, 0);
            currentDisplacements.Add(new Vector(incrementalNodeDisplacements));
            new Vector(incrementalNodeDisplacements).CopyTo(displacementsOfCurrentIncrement.Data, 0);

            double currentDisplacementX_A = currentDisplacements[0];
            double currentlDisplacementY_A = currentDisplacements[1];
            double currentRotationZ_A = currentDisplacements[2];
            double currentDisplacementX_B = currentDisplacements[3];
            double currentDisplacementY_B = currentDisplacements[4];
            double currentRotationZ_Β = currentDisplacements[5];

            double dX = ((nodes[1].X - nodes[0].X) + this.currentDisplacements[0]) - this.currentDisplacements[3];
            double dY = ((nodes[1].Y - nodes[0].Y) + this.currentDisplacements[1]) - this.currentDisplacements[4];
            this.currentLength = Math.Sqrt((dX * dX) + (dY * dY));
            double axisAngle = 0d;
            if ((dY == 0) && (dX == currentLength)) { axisAngle = 0d; }
            else
                if ((dY == 0) && (dX == -currentLength)) { axisAngle = Math.PI; }
            else { axisAngle = 2.0 * Math.Atan((currentLength - dX) / dY); }

            double initialdX = nodes[1].X - nodes[0].X;
            double initialdY = nodes[1].Y - nodes[0].Y;
            double initialLength = Math.Sqrt((initialdX * initialdX) + (initialdY * initialdY));
            double axisAngleInitial = 0d;
            if ((initialdY == 0) && (initialdX == initialLength)) { axisAngleInitial = 0d; }
            else
                if ((initialdY == 0) && (initialdX == -initialLength)) { axisAngleInitial = Math.PI; }
            else { axisAngleInitial = 2.0 * Math.Atan((initialLength - dX) / dY); }
            
            double symmetricAngle = currentRotationZ_Β - currentRotationZ_A;
            double antiSymmetricAngle = currentRotationZ_Β + currentRotationZ_A - 2.0 * (axisAngle - axisAngleInitial);

            currentRotationMatrix = RotationMatrix.CalculateRotationMatrixBeam2D(axisAngle);
            double extension = this.currentLength - this.initialLength;
            this.naturalDeformations[0] = extension;
            this.naturalDeformations[1] = symmetricAngle;
            this.naturalDeformations[2] = antiSymmetricAngle;
        }

        private void CalculateUpdatedBeamAxis()
        {
            //currentRotationMatrix.Multiply(initialAxisX, beamAxisX.Data);
            //currentRotationMatrix.Multiply(initialAxisY, beamAxisY.Data);
            //currentRotationMatrix.Multiply(initialAxisZ, beamAxisZ.Data);

            //double dX = ((nodes[1].X - nodes[0].X) + this.currentDisplacements[6]) - this.currentDisplacements[0];
            //double dY = ((nodes[1].Y - nodes[0].Y) + this.currentDisplacements[7]) - this.currentDisplacements[1];
            //double dZ = ((nodes[1].Z - nodes[0].Z) + this.currentDisplacements[8]) - this.currentDisplacements[2];
            //this.currentLength = Math.Sqrt((dX * dX) + (dY * dY) + (dZ * dZ));
            //var delta = new Vector(new[] { dX, dY, dZ });
            //this.currentBeamAxis[0] = dX / currentLength + beamAxisX[0];
            //this.currentBeamAxis[1] = dY / currentLength + beamAxisX[1];
            //this.currentBeamAxis[2] = dZ / currentLength + beamAxisX[2];
            //this.currentBeamAxis.Scale(1d / this.currentBeamAxis.Norm);
        }

        private void InitializeElementAxes()
        {
            double deltaX = nodes[1].X - nodes[0].X;
            double deltaY = nodes[1].Y - nodes[0].Y;
            this.currentLength = Math.Sqrt((deltaX * deltaX) + (deltaY * deltaY));
            double axisAngle = 0d;
            if ((deltaY == 0) && (deltaX == currentLength)) { axisAngle = 0d; }
            else
                if ((deltaY == 0) && (deltaX == -currentLength)) { axisAngle = Math.PI; }
            else { axisAngle = 2.0 * Math.Atan((currentLength - deltaX) / deltaY); }
            
            currentRotationMatrix = RotationMatrix.CalculateRotationMatrixBeam2D(axisAngle);
        }

        private void UpdateNaturalDeformations(Vector vectorPartDifference)
        {
            double extension = this.currentLength - this.initialLength;
            Matrix2D currentRotationMatrixTranspose = currentRotationMatrix.Transpose();
            var symmetricRotation = new Vector(AXIS_COUNT);

            currentRotationMatrixTranspose.Multiply(vectorPartDifference, symmetricRotation.Data);

            var antisymmetricRotation = new Vector(AXIS_COUNT);

            currentRotationMatrixTranspose.Multiply(this.beamAxisX ^ this.currentBeamAxis, antisymmetricRotation.Data);
            this.naturalDeformations[0] = symmetricRotation[0] * 4.0;
            this.naturalDeformations[1] = symmetricRotation[1] * 4.0;
            this.naturalDeformations[2] = symmetricRotation[2] * 4.0;
            this.naturalDeformations[3] = extension;
            this.naturalDeformations[4] = antisymmetricRotation[1] * 4.0;
            this.naturalDeformations[5] = antisymmetricRotation[2] * 4.0;
        }

        private void UpdateRotationMatrix()
        {
            //currentRotationMatrix = RotationMatrix.CalculateRotationMatrixBeam2D(this.axisAngle);
        }
    }
}
