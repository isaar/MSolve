using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.IGA.Elements;
using ISAAR.MSolve.IGA.Problems.Structural.Elements;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Materials.Interfaces;

namespace ISAAR.MSolve.IGA.Entities
{
    public class CollocationPatch:IAsymmetricSubdomain
    {
        private readonly Dictionary<int, Edge> edgesDictionary = new Dictionary<int, Edge>();
        private readonly Dictionary<int, Face> facesDictionary = new Dictionary<int, Face>();

        private readonly List<ControlPoint> controlPoints = new List<ControlPoint>();

        public Table<INode, IDofType, double> Constraints { get; } = new Table<INode, IDofType, double>();

        IReadOnlyList<IElement> ISubdomain.Elements => Elements;
        public List<ICollocationElement> Elements { get; } = new List<ICollocationElement>();

        public int ID { get; }

        IReadOnlyList<INode> ISubdomain.Nodes => controlPoints;
        public IReadOnlyList<ControlPoint> ControlPoints => controlPoints;

        public ISubdomainConstrainedDofOrdering ConstrainedDofOrdering { get; set; }
        public ISubdomainFreeDofOrdering FreeDofOrdering { get; set; }

        public Vector Forces { get; set; }

        public bool StiffnessModified { get; set; }
        public bool ConnectivityModified { get; set; }

        public Table<ControlPoint, IDofType, double> ControlPointLoads { get; set; }

        public Dictionary<int, Edge> EdgesDictionary
        {
            get { return edgesDictionary; }
        }

        public Dictionary<int, Face> FacesDictionary
        {
            get { return facesDictionary; }
        }

        public double[] CalculateElementIncrementalConstraintDisplacements(IElement element, double constraintScalingFactor)
        {
            var elementNodalDisplacements = new double[FreeDofOrdering.CountElementDofs(element)];
            ApplyConstraintDisplacements(element, elementNodalDisplacements, Constraints);
            return elementNodalDisplacements;
        }

        private static void ApplyConstraintDisplacements(IElement element, double[] elementNodalDisplacements,
            Table<INode, IDofType, double> constraints)
        {
            int elementDofIdx = 0;
            IReadOnlyList<INode> nodes = element.ElementType.DofEnumerator.GetNodesForMatrixAssembly(element);
            IReadOnlyList<IReadOnlyList<IDofType>> dofs = element.ElementType.DofEnumerator.GetDofTypesForMatrixAssembly(element);
            for (int i = 0; i < nodes.Count; ++i)
            {
                //bool isConstrainedNode = constraintsDictionary.TryGetValue(nodes[i].ID, 
                //    out Dictionary<DOFType, double> constrainedDOFs);
                bool isConstrainedNode = constraints.TryGetDataOfRow(nodes[i],
                    out IReadOnlyDictionary<IDofType, double> constrainedDOFs);
                if (isConstrainedNode)
                {
                    foreach (IDofType dofType in dofs[i])
                    {
                        bool isConstrainedDof = constrainedDOFs.TryGetValue(dofType, out double constraintDisplacement);
                        //if (isConstrainedNode && isConstrainedDof)
                        if (isConstrainedDof)
                        {
                            Debug.Assert(elementNodalDisplacements[elementDofIdx] == 0); // TODO: and why is this an assumption?
                            elementNodalDisplacements[elementDofIdx] = constraintDisplacement;
                        }
                        ++elementDofIdx;
                    }
                }
                else elementDofIdx += dofs[i].Count;
            }
        }

        public double[] CalculateElementDisplacements(Element element, IVectorView globalDisplacementVector)//QUESTION: would it be maybe more clear if we passed the constraintsDictionary as argument??
        {
            var elementNodalDisplacements = new double[FreeDofOrdering.CountElementDofs(element)];
            FreeDofOrdering.ExtractVectorElementFromSubdomain(element, globalDisplacementVector);
            ApplyConstraintDisplacements(element, elementNodalDisplacements, Constraints);
            return elementNodalDisplacements;
        }

        public void ClearMaterialStresses()
        {
            //foreach (Element element in Elements) element.ElementType.ClearMaterialStresses();
        }

        public void DefineControlPointsFromElements()
        {
            var cpComparer = Comparer<ControlPoint>.Create((node1, node2) => node1.ID - node2.ID);
            var cpSet = new SortedSet<ControlPoint>(cpComparer);
            foreach (Element element in Elements)
            {
                foreach (ControlPoint node in element.ControlPoints) cpSet.Add(node);
            }
            controlPoints.AddRange(cpSet);
        }

        public void ExtractConstraintsFromGlobal(Table<INode, IDofType, double> globalConstraints)
        {
            foreach (ControlPoint controlPoint in ControlPoints)
            {
                bool isControlPointConstrained = globalConstraints.TryGetDataOfRow(controlPoint,
                    out IReadOnlyDictionary<IDofType, double> constraintsOfNode);
                if (isControlPointConstrained)
                {
                    foreach (var dofDisplacementPair in constraintsOfNode)
                    {
                        Constraints[controlPoint, dofDisplacementPair.Key] = dofDisplacementPair.Value;
                    }
                }
            }

            // This is probably faster but assumes that nodes store their prescribed displacements, which I hate.
            //foreach (Node node in Nodes)
            //{
            //    if (node.Constraints == null) continue;
            //    foreach (Constraint constraint in node.Constraints) Constraints[node, constraint.DOF] = constraint.Amount;
            //}
        }

        public IVector GetRhsFromSolution(IVectorView solution, IVectorView dSolution)
        {
            var forces = Vector.CreateZero(FreeDofOrdering.NumFreeDofs); //TODO: use Vector
            foreach (Element element in Elements)
            {
                double[] localSolution = CalculateElementDisplacements(element, solution);
                double[] localdSolution = CalculateElementDisplacements(element, dSolution);
                element.ElementType.CalculateStresses(element, localSolution, localdSolution);
                if (element.ElementType.MaterialModified)
                    element.Patch.StiffnessModified = true;
                var f = element.ElementType.CalculateForces(element, localSolution, localdSolution);
                FreeDofOrdering.AddVectorElementToSubdomain(element, f, forces);
            }
            return forces;
        }

        public void ResetMaterialsModifiedProperty()
        {
            this.StiffnessModified = false;
            foreach (Element element in Elements) element.ElementType.ResetMaterialModified();

        }

        public void ScaleConstraints(double scalingFactor) => Constraints.ModifyValues((u) => scalingFactor * u);

        public void SaveMaterialState()
        {
            //foreach (Element element in Elements) element.ElementType.SaveMaterialState();
        }

        #region PatchData
        public int NumberOfDimensions { get; set; }
        public double Thickness { get; set; }
        public IFiniteElementMaterial Material { get; set; }

        public int NumberOfControlPointsKsi { get; set; }
        public int NumberOfControlPointsHeta { get; set; }
        public int NumberOfControlPointsZeta { get; set; }

        public int DegreeKsi { get; set; }
        public int DegreeHeta { get; set; }
        public int DegreeZeta { get; set; }

        public Vector KnotValueVectorKsi { get; set; }
        public Vector KnotValueVectorHeta { get; set; }
        public Vector KnotValueVectorZeta { get; set; }
        public ISubdomainFreeDofOrdering FreeDofRowOrdering { get; set; }
        public ISubdomainFreeDofOrdering FreeDofColOrdering { get; set; }

        public ISubdomainConstrainedDofOrdering ConstrainedDofRowOrdering { get; set; }
        public ISubdomainConstrainedDofOrdering ConstrainedDofColOrdering { get; set; }
        #endregion

        public void CreateCollocationPatchData()
        {
            if (this.NumberOfDimensions == 2)
            {
                CreatePatchCollocationData2D();
            }
            else
            {
                CreatePatchCollocationData3D();
            }
        }

        private void CreatePatchCollocationData3D()
        {
            CreateCollocationPoints3D();
        }

        private void CreateCollocationPoints3D()
        {
            #region Collocation Points
            var collocationPoints = new List<CollocationPoint3D>();
            var index = 0;
            for (int i = 0; i < NumberOfControlPointsKsi; i++)
            {
                var coordinateKsi = 0.0;
                for (int j = 1; j <= DegreeKsi; j++)
                    coordinateKsi += KnotValueVectorKsi[i + j];
                coordinateKsi /= DegreeKsi;

                for (int j = 0; j < NumberOfControlPointsHeta; j++)
                {
                    var coordinateHeta = 0.0;
                    for (int k = 1; k <= DegreeHeta; k++)
                        coordinateHeta += KnotValueVectorHeta[j + k];

                    coordinateHeta /= DegreeHeta;

                    for (int k = 0; k < NumberOfControlPointsZeta; k++)
                    {
                        var coordinateZeta = 0.0;
                        for (int m = 1; m <= DegreeZeta; m++)
                            coordinateZeta += KnotValueVectorZeta[k + m];

                        coordinateZeta /= DegreeZeta;

                        var isBoundary = (i == 0 || i == NumberOfControlPointsKsi - 1 || j == 0 ||
                                          j == NumberOfControlPointsHeta - 1 || k == 0 || k == NumberOfControlPointsZeta - 1);
                        var collocationPoint = new CollocationPoint3D(index++, coordinateKsi, coordinateHeta,
                            coordinateZeta, isBoundary);
                        if (i == 0) collocationPoint.Surfaces.Add(Surface.Left);
                        if (i == NumberOfControlPointsKsi - 1) collocationPoint.Surfaces.Add(Surface.Right);
                        if (j == 0) collocationPoint.Surfaces.Add(Surface.Front);
                        if (j == NumberOfControlPointsHeta - 1) collocationPoint.Surfaces.Add(Surface.Back);
                        if (k == 0) collocationPoint.Surfaces.Add(Surface.Bottom);
                        if (k == NumberOfControlPointsZeta - 1) collocationPoint.Surfaces.Add(Surface.Top);
                        collocationPoints.Add(collocationPoint);
                    }
                }
            }
            #endregion

            #region Knots
            Vector singleKnotValuesKsi = KnotValueVectorKsi.RemoveDuplicatesFindMultiplicity()[0];
            Vector singleKnotValuesHeta = KnotValueVectorHeta.RemoveDuplicatesFindMultiplicity()[0];
            Vector singleKnotValuesZeta = KnotValueVectorZeta.RemoveDuplicatesFindMultiplicity()[0];

            List<Knot> knots = new List<Knot>();

            int id = 0;
            for (int i = 0; i < singleKnotValuesKsi.Length; i++)
            {
                for (int j = 0; j < singleKnotValuesHeta.Length; j++)
                {
                    for (int k = 0; k < singleKnotValuesZeta.Length; k++)
                    {
                        knots.Add(new Knot() { ID = id, Ksi = singleKnotValuesKsi[i], Heta = singleKnotValuesHeta[j], Zeta = singleKnotValuesZeta[k] });
                        id++;
                    }

                }
            }
            #endregion

            #region Elements
            Vector singlesKnotValuesKsi = KnotValueVectorKsi.RemoveDuplicatesFindMultiplicity()[0];
            Vector multiplicityKsi = KnotValueVectorKsi.RemoveDuplicatesFindMultiplicity()[1];
            Vector singlesKnotValuesHeta = KnotValueVectorHeta.RemoveDuplicatesFindMultiplicity()[0];
            Vector multiplicityHeta = KnotValueVectorHeta.RemoveDuplicatesFindMultiplicity()[1];
            Vector singlesKnotValuesZeta = KnotValueVectorZeta.RemoveDuplicatesFindMultiplicity()[0];
            Vector multiplicityZeta = KnotValueVectorZeta.RemoveDuplicatesFindMultiplicity()[1];

            int numberOfElementsKsi = singlesKnotValuesKsi.Length - 1;
            int numberOfElementsHeta = singlesKnotValuesHeta.Length - 1;
            int numberOfElementsZeta = singlesKnotValuesZeta.Length - 1;

            if (numberOfElementsKsi * numberOfElementsHeta * numberOfElementsZeta == 0)
            {
                throw new NullReferenceException("Number of Elements should be defined before Element Connectivity");
            }

            for (int i = 0; i < numberOfElementsKsi; i++)
            {
                for (int j = 0; j < numberOfElementsHeta; j++)
                {
                    for (int k = 0; k < numberOfElementsZeta; k++)
                    {
                        IList<Knot> knotsOfElement = new List<Knot>();
                        knotsOfElement.Add(knots[i * singlesKnotValuesZeta.Length * singlesKnotValuesHeta.Length + j * singlesKnotValuesZeta.Length + k]);
                        knotsOfElement.Add(knots[i * singlesKnotValuesZeta.Length * singlesKnotValuesHeta.Length + j * singlesKnotValuesZeta.Length + k + 1]);
                        knotsOfElement.Add(knots[i * singlesKnotValuesZeta.Length * singlesKnotValuesHeta.Length + (j + 1) * singlesKnotValuesZeta.Length + k]);
                        knotsOfElement.Add(knots[i * singlesKnotValuesZeta.Length * singlesKnotValuesHeta.Length + (j + 1) * singlesKnotValuesZeta.Length + k + 1]);
                        knotsOfElement.Add(knots[(i + 1) * singlesKnotValuesZeta.Length * singlesKnotValuesHeta.Length + j * singlesKnotValuesZeta.Length + k]);
                        knotsOfElement.Add(knots[(i + 1) * singlesKnotValuesZeta.Length * singlesKnotValuesHeta.Length + j * singlesKnotValuesZeta.Length + k + 1]);
                        knotsOfElement.Add(knots[(i + 1) * singlesKnotValuesZeta.Length * singlesKnotValuesHeta.Length + (j + 1) * singlesKnotValuesZeta.Length + k]);
                        knotsOfElement.Add(knots[(i + 1) * singlesKnotValuesZeta.Length * singlesKnotValuesHeta.Length + (j + 1) * singlesKnotValuesZeta.Length + k + 1]);

                        int multiplicityElementKsi = 0;
                        if (multiplicityKsi[i + 1] - this.DegreeKsi > 0)
                        {
                            multiplicityElementKsi = (int)multiplicityKsi[i + 1] - DegreeKsi;
                        }

                        int multiplicityElementHeta = 0;
                        if (multiplicityHeta[j + 1] - this.DegreeHeta > 0)
                        {
                            multiplicityElementHeta = (int)multiplicityHeta[j + 1] - this.DegreeHeta;
                        }

                        int multiplicityElementZeta = 0;
                        if (multiplicityZeta[k + 1] - this.DegreeZeta > 0)
                        {
                            multiplicityElementZeta = (int)multiplicityZeta[k + 1] - this.DegreeZeta;
                        }

                        int nurbsSupportKsi = this.DegreeKsi + 1;
                        int nurbsSupportHeta = this.DegreeHeta + 1;
                        int nurbsSupportZeta = this.DegreeZeta + 1;

                        IList<ControlPoint> elementControlPoints = new List<ControlPoint>();

                        for (int l = 0; l < nurbsSupportKsi; l++)
                        {
                            for (int m = 0; m < nurbsSupportHeta; m++)
                            {
                                for (int n = 0; n < nurbsSupportZeta; n++)
                                {
                                    int controlPointID = (i + multiplicityElementKsi) * NumberOfControlPointsHeta * NumberOfControlPointsZeta +
                                        (j + multiplicityElementHeta) * NumberOfControlPointsZeta + (k + multiplicityElementZeta) +
                                        l * NumberOfControlPointsHeta * NumberOfControlPointsZeta + m * NumberOfControlPointsZeta + n;

                                    elementControlPoints.Add(ControlPoints[controlPointID]);
                                }
                            }
                        }

                        var elementCollocationPoints = collocationPoints.Where(cp =>
                            knotsOfElement[0].Ksi <= cp.Xi && cp.Xi <= knotsOfElement[7].Ksi &&
                            knotsOfElement[0].Heta <= cp.Eta && cp.Eta <= knotsOfElement[7].Heta &&
                            knotsOfElement[0].Zeta <= cp.Zeta && cp.Zeta <= knotsOfElement[7].Zeta).ToList();

                        foreach (var point in elementCollocationPoints)
                        {
                            if (Elements.Any(e => e.ID == point.ID)) continue;
                            var element = new NURBSElement3DCollocation
                            {
                                ID = point.ID,
                                Patch = this,
                                ElementType = new NURBSElement3DCollocation(),
                                CollocationPoint = point
                            };
                            element.AddKnots(knotsOfElement);
                            element.AddControlPoints(elementControlPoints.ToList<ControlPoint>());
                            Elements.Add(element);
                        }
                    }
                }
            }

            var orderedElements = Elements.OrderBy(e => ((ICollocationElement)e).CollocationPoint.ID).ToList();
            Elements.Clear();
            Elements.AddRange(orderedElements);

            #endregion

        }

        private void CreatePatchCollocationData2D()
        {
            CreateCollocationPoints2D();
        }

        private void CreateCollocationPoints2D()
        {
            #region CollocationPoints
            var collocationPoints = new List<CollocationPoint2D>();
            var index = 0;
            for (int i = 0; i < NumberOfControlPointsKsi; i++)
            {
                var coordinateKsi = 0.0;
                for (int j = 1; j <= DegreeKsi; j++)
                    coordinateKsi += KnotValueVectorKsi[i + j];
                coordinateKsi /= DegreeKsi;

                for (int j = 0; j < NumberOfControlPointsHeta; j++)
                {
                    var coordinateHeta = 0.0;
                    for (int k = 1; k <= DegreeHeta; k++)
                        coordinateHeta += KnotValueVectorHeta[j + k];

                    coordinateHeta /= DegreeHeta;

                    var isBoundary = (i == 0 || i == NumberOfControlPointsKsi - 1 || j == 0 ||
                                      j == NumberOfControlPointsHeta - 1);
                    collocationPoints.Add(new CollocationPoint2D(index++, coordinateKsi, coordinateHeta, isBoundary));
                }
            }
            #endregion

            #region Knots
            Vector singleKnotValuesKsi = KnotValueVectorKsi.RemoveDuplicatesFindMultiplicity()[0];
            Vector singleKnotValuesHeta = KnotValueVectorHeta.RemoveDuplicatesFindMultiplicity()[0];

            List<Knot> knots = new List<Knot>();

            int id = 0;
            for (int i = 0; i < singleKnotValuesKsi.Length; i++)
            {
                for (int j = 0; j < singleKnotValuesHeta.Length; j++)
                {
                    knots.Add(new Knot() { ID = id, Ksi = singleKnotValuesKsi[i], Heta = singleKnotValuesHeta[j], Zeta = 0.0 });
                    id++;
                }
            }
            #endregion

            #region Elements
            Vector multiplicityKsi = KnotValueVectorKsi.RemoveDuplicatesFindMultiplicity()[1];
            Vector multiplicityHeta = KnotValueVectorHeta.RemoveDuplicatesFindMultiplicity()[1];

            int numberOfElementsKsi = singleKnotValuesKsi.Length - 1;
            int numberOfElementsHeta = singleKnotValuesHeta.Length - 1;
            if (numberOfElementsKsi * numberOfElementsHeta == 0)
            {
                throw new NullReferenceException("Number of Elements should be defined before Element Connectivity");
            }

            for (int i = 0; i < numberOfElementsKsi; i++)
            {
                for (int j = 0; j < numberOfElementsHeta; j++)
                {
                    IList<Knot> knotsOfElement = new List<Knot>();
                    knotsOfElement.Add(knots[i * singleKnotValuesHeta.Length + j]);
                    knotsOfElement.Add(knots[i * singleKnotValuesHeta.Length + j + 1]);
                    knotsOfElement.Add(knots[(i + 1) * singleKnotValuesHeta.Length + j]);
                    knotsOfElement.Add(knots[(i + 1) * singleKnotValuesHeta.Length + j + 1]);

                    int multiplicityElementKsi = 0;
                    if (multiplicityKsi[i + 1] - this.DegreeKsi > 0)
                    {
                        multiplicityElementKsi = (int)multiplicityKsi[i + 1] - this.DegreeKsi;
                    }

                    int multiplicityElementHeta = 0;
                    if (multiplicityHeta[j + 1] - this.DegreeHeta > 0)
                    {
                        multiplicityElementHeta = (int)multiplicityHeta[j + 1] - this.DegreeHeta;
                    }

                    int nurbsSupportKsi = this.DegreeKsi + 1;
                    int nurbsSupportHeta = this.DegreeHeta + 1;

                    IList<ControlPoint> elementControlPoints = new List<ControlPoint>();

                    for (int k = 0; k < nurbsSupportKsi; k++)
                    {
                        for (int l = 0; l < nurbsSupportHeta; l++)
                        {
                            int controlPointID = (i + multiplicityElementKsi) * this.NumberOfControlPointsHeta +
                                (j + multiplicityElementHeta) + k * this.NumberOfControlPointsHeta + l;
                            elementControlPoints.Add(this.ControlPoints[controlPointID]);
                        }
                    }

                    var elementCollocationPoints = collocationPoints.Where(cp =>
                        knotsOfElement[0].Ksi <= cp.Xi && cp.Xi <= knotsOfElement[3].Ksi &&
                        knotsOfElement[0].Heta <= cp.Eta && cp.Eta <= knotsOfElement[3].Heta).ToList();

                    foreach (var point in elementCollocationPoints)
                    {
                        if (Elements.Any(e => e.ID == point.ID)) continue;
                        var element = new NURBSElement2DCollocation
                        {
                            ID = point.ID,
                            Patch = this,
                            ElementType = new NURBSElement2DCollocation(),
                            CollocationPoint = point
                        };
                        element.AddKnots(knotsOfElement);
                        element.AddControlPoints(elementControlPoints.ToList<ControlPoint>());
                        Elements.Add(element);
                    }
                }
            }

            var orderedElements = Elements.OrderBy(e => ((ICollocationElement)e).CollocationPoint.ID).ToList();
            Elements.Clear();
            Elements.AddRange(orderedElements);


            #endregion
        }

    }

}

