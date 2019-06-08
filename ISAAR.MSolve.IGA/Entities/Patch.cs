using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using ISAAR.MSolve.Discretization.Commons;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Geometry.Coordinates;
using ISAAR.MSolve.IGA.Elements;
using ISAAR.MSolve.IGA.Problems.Structural.Elements;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.Materials.Interfaces;

namespace ISAAR.MSolve.IGA.Entities
{
	public class Patch: ISubdomain
    {
		private readonly Dictionary<int, Edge> edgesDictionary = new Dictionary<int, Edge>();
		private readonly Dictionary<int, Face> facesDictionary = new Dictionary<int, Face>();
		
		private readonly List<ControlPoint> controlPoints = new List<ControlPoint>();

		public Table<INode, IDofType, double> Constraints { get; } = new Table<INode, IDofType, double>();

		IReadOnlyList<IElement> ISubdomain.Elements => Elements;
		public List<Element> Elements { get; } = new List<Element>();
		
		public int ID { get; }

		IReadOnlyList<INode> ISubdomain.Nodes => controlPoints;
		public IReadOnlyList<ControlPoint> ControlPoints => controlPoints;

        public ISubdomainConstrainedDofOrdering ConstrainedDofOrdering { get; set; }
        public ISubdomainFreeDofOrdering FreeDofOrdering { get; set; }

		public Vector Forces { get; set; }

		public bool StiffnessModified { get; set; }
        public bool ConnectivityModified { get; set; } = true; // At first it is modified

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
            SubdomainConstrainedDofOrderingBase.ApplyConstraintDisplacements(element, elementNodalDisplacements, Constraints);
			return elementNodalDisplacements;
		}

		public double[] CalculateElementDisplacements(Element element, IVectorView globalDisplacementVector)//QUESTION: would it be maybe more clear if we passed the constraintsDictionary as argument??
		{
			var elementNodalDisplacements = new double[FreeDofOrdering.CountElementDofs(element)];
			FreeDofOrdering.ExtractVectorElementFromSubdomain(element, globalDisplacementVector);
            SubdomainConstrainedDofOrderingBase.ApplyConstraintDisplacements(element, elementNodalDisplacements, Constraints);
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


        public void CreatePatchData()
		{
			if (this.NumberOfDimensions == 2)
			{
				CreatePatchData2D();
			}
			else
			{
				CreatePatchData3D();
			}
		}

//		public void CreateCollocationPatchData()
//		{
//			if (this.NumberOfDimensions == 2)
//			{
//				CreatePatchCollocationData2D();
//			}
//			else
//			{
//				CreatePatchCollocationData3D();
//			}
//		}

//        private void CreatePatchCollocationData3D()
//        {
//            CreateCollocationPoints3D();
//        }

//        private void CreateCollocationPoints3D()
//        {
//            #region Collocation Points
//            var collocationPoints = new List<CollocationPoint3D>();
//            var index = 0;
//            for (int i = 0; i < NumberOfControlPointsKsi; i++)
//            {
//                var coordinateKsi = 0.0;
//                for (int j = 1; j <= DegreeKsi; j++)
//                    coordinateKsi += KnotValueVectorKsi[i + j];
//                coordinateKsi /= DegreeKsi;

//                for (int j = 0; j < NumberOfControlPointsHeta; j++)
//                {
//                    var coordinateHeta = 0.0;
//                    for (int k = 1; k <= DegreeHeta; k++)
//                        coordinateHeta += KnotValueVectorHeta[j + k];

//                    coordinateHeta /= DegreeHeta;

//                    for (int k = 0; k < NumberOfControlPointsZeta; k++)
//                    {
//                        var coordinateZeta = 0.0;
//                        for (int m = 1; m <= DegreeZeta; m++)
//                            coordinateZeta += KnotValueVectorZeta[k + m];

//                        coordinateZeta /= DegreeZeta;

//                        var isBoundary = (i == 0 || i == NumberOfControlPointsKsi - 1 || j == 0 ||
//                                          j == NumberOfControlPointsHeta - 1 || k == 0 || k == NumberOfControlPointsZeta-1);
//                        var collocationPoint = new CollocationPoint3D(index++, coordinateKsi, coordinateHeta,
//                            coordinateZeta, isBoundary);
//                        if (i == 0) collocationPoint.Surfaces.Add(Surface.Left);
//                        if (i == NumberOfControlPointsKsi - 1) collocationPoint.Surfaces.Add(Surface.Right);
//                        if (j == 0) collocationPoint.Surfaces.Add(Surface.Front);
//                        if (j == NumberOfControlPointsHeta - 1) collocationPoint.Surfaces.Add(Surface.Back);
//                        if (k == 0) collocationPoint.Surfaces.Add(Surface.Bottom);
//                        if (k == NumberOfControlPointsZeta - 1) collocationPoint.Surfaces.Add(Surface.Top);
//                        collocationPoints.Add(collocationPoint);
//                    }
//                }
//            }
//            #endregion
            
//            #region Knots
//            Vector singleKnotValuesKsi = KnotValueVectorKsi.RemoveDuplicatesFindMultiplicity()[0];
//            Vector singleKnotValuesHeta = KnotValueVectorHeta.RemoveDuplicatesFindMultiplicity()[0];
//            Vector singleKnotValuesZeta = KnotValueVectorZeta.RemoveDuplicatesFindMultiplicity()[0];

//            List<Knot> knots = new List<Knot>();

//            int id = 0;
//            for (int i = 0; i < singleKnotValuesKsi.Length; i++)
//            {
//                for (int j = 0; j < singleKnotValuesHeta.Length; j++)
//                {
//                    for (int k = 0; k < singleKnotValuesZeta.Length; k++)
//                    {
//                        knots.Add(new Knot() { ID = id, Ksi = singleKnotValuesKsi[i], Heta = singleKnotValuesHeta[j], Zeta = singleKnotValuesZeta[k] });
//                        id++;
//                    }

//                }
//            }
//            #endregion

//            #region Elements
//            Vector singlesKnotValuesKsi = KnotValueVectorKsi.RemoveDuplicatesFindMultiplicity()[0];
//            Vector multiplicityKsi = KnotValueVectorKsi.RemoveDuplicatesFindMultiplicity()[1];
//            Vector singlesKnotValuesHeta = KnotValueVectorHeta.RemoveDuplicatesFindMultiplicity()[0];
//            Vector multiplicityHeta = KnotValueVectorHeta.RemoveDuplicatesFindMultiplicity()[1];
//            Vector singlesKnotValuesZeta = KnotValueVectorZeta.RemoveDuplicatesFindMultiplicity()[0];
//            Vector multiplicityZeta = KnotValueVectorZeta.RemoveDuplicatesFindMultiplicity()[1];

//            int numberOfElementsKsi = singlesKnotValuesKsi.Length - 1;
//            int numberOfElementsHeta = singlesKnotValuesHeta.Length - 1;
//            int numberOfElementsZeta = singlesKnotValuesZeta.Length - 1;

//            if (numberOfElementsKsi * numberOfElementsHeta * numberOfElementsZeta == 0)
//            {
//                throw new NullReferenceException("Number of Elements should be defined before Element Connectivity");
//            }

//            for (int i = 0; i < numberOfElementsKsi; i++)
//            {
//                for (int j = 0; j < numberOfElementsHeta; j++)
//                {
//                    for (int k = 0; k < numberOfElementsZeta; k++)
//                    {
//                        IList<Knot> knotsOfElement = new List<Knot>();
//                        knotsOfElement.Add(knots[i * singlesKnotValuesZeta.Length * singlesKnotValuesHeta.Length + j * singlesKnotValuesZeta.Length + k]);
//                        knotsOfElement.Add(knots[i * singlesKnotValuesZeta.Length * singlesKnotValuesHeta.Length + j * singlesKnotValuesZeta.Length + k + 1]);
//                        knotsOfElement.Add(knots[i * singlesKnotValuesZeta.Length * singlesKnotValuesHeta.Length + (j + 1) * singlesKnotValuesZeta.Length + k]);
//                        knotsOfElement.Add(knots[i * singlesKnotValuesZeta.Length * singlesKnotValuesHeta.Length + (j + 1) * singlesKnotValuesZeta.Length + k + 1]);
//                        knotsOfElement.Add(knots[(i + 1) * singlesKnotValuesZeta.Length * singlesKnotValuesHeta.Length + j * singlesKnotValuesZeta.Length + k]);
//                        knotsOfElement.Add(knots[(i + 1) * singlesKnotValuesZeta.Length * singlesKnotValuesHeta.Length + j * singlesKnotValuesZeta.Length + k + 1]);
//                        knotsOfElement.Add(knots[(i + 1) * singlesKnotValuesZeta.Length * singlesKnotValuesHeta.Length + (j + 1) * singlesKnotValuesZeta.Length + k]);
//                        knotsOfElement.Add(knots[(i + 1) * singlesKnotValuesZeta.Length * singlesKnotValuesHeta.Length + (j + 1) * singlesKnotValuesZeta.Length + k + 1]);

//                        int multiplicityElementKsi = 0;
//                        if (multiplicityKsi[i + 1] - this.DegreeKsi > 0)
//                        {
//                            multiplicityElementKsi = (int)multiplicityKsi[i + 1] - DegreeKsi;
//                        }

//                        int multiplicityElementHeta = 0;
//                        if (multiplicityHeta[j + 1] - this.DegreeHeta > 0)
//                        {
//                            multiplicityElementHeta = (int)multiplicityHeta[j + 1] - this.DegreeHeta;
//                        }

//                        int multiplicityElementZeta = 0;
//                        if (multiplicityZeta[k + 1] - this.DegreeZeta > 0)
//                        {
//                            multiplicityElementZeta = (int)multiplicityZeta[k + 1] - this.DegreeZeta;
//                        }

//                        int nurbsSupportKsi = this.DegreeKsi + 1;
//                        int nurbsSupportHeta = this.DegreeHeta + 1;
//                        int nurbsSupportZeta = this.DegreeZeta + 1;

//                        IList<ControlPoint> elementControlPoints = new List<ControlPoint>();

//                        for (int l = 0; l < nurbsSupportKsi; l++)
//                        {
//                            for (int m = 0; m < nurbsSupportHeta; m++)
//                            {
//                                for (int n = 0; n < nurbsSupportZeta; n++)
//                                {
//                                    int controlPointID = (i + multiplicityElementKsi) * NumberOfControlPointsHeta * NumberOfControlPointsZeta +
//                                        (j + multiplicityElementHeta) * NumberOfControlPointsZeta + (k + multiplicityElementZeta) +
//                                        l * NumberOfControlPointsHeta * NumberOfControlPointsZeta + m * NumberOfControlPointsZeta + n;

//                                    elementControlPoints.Add(ControlPoints[controlPointID]);
//                                }
//                            }
//                        }

//                        var elementCollocationPoints = collocationPoints.Where(cp =>
//                            knotsOfElement[0].Ksi <= cp.Xi && cp.Xi <= knotsOfElement[7].Ksi &&
//                            knotsOfElement[0].Heta <= cp.Eta && cp.Eta <= knotsOfElement[7].Heta &&
//                            knotsOfElement[0].Zeta <= cp.Zeta && cp.Zeta <= knotsOfElement[7].Zeta).ToList();

//                        foreach (var point in elementCollocationPoints)
//                        {
//                            if (Elements.Any(e => e.ID == point.ID)) continue;
//                            Element element = new NURBSElement3DCollocation
//                            {
//                                ID = point.ID,
//                                Patch = this,
//                                ElementType = new NURBSElement3DCollocation(),
//                                CollocationPoint = point
//                            };
//                            element.AddKnots(knotsOfElement);
//                            element.AddControlPoints(elementControlPoints.ToList<ControlPoint>());
//                            Elements.Add(element);
//                        }
//                    }
//                }
//            }

//            var orderedElements = Elements.OrderBy(e => ((ICollocationElement) e).CollocationPoint.ID).ToList();
//            Elements.Clear();
//            Elements.AddRange(orderedElements);

//            #endregion

//        }

//        private void CreatePatchCollocationData2D()
//		{
//			CreateCollocationPoints2D();
//		}

//		private void CreateCollocationPoints2D()
//		{
//			#region CollocationPoints
//			var collocationPoints= new List<CollocationPoint2D>();
//			var index = 0;
//			for (int i = 0; i < NumberOfControlPointsKsi; i++)
//			{
//				var coordinateKsi = 0.0;
//				for (int j = 1; j <= DegreeKsi; j++)
//					coordinateKsi += KnotValueVectorKsi[i+j];
//				coordinateKsi /= DegreeKsi;

//				for (int j = 0; j < NumberOfControlPointsHeta; j++)
//				{
//					var coordinateHeta = 0.0;
//					for (int k = 1; k <= DegreeHeta; k++)
//						coordinateHeta += KnotValueVectorHeta[j + k];

//					coordinateHeta /= DegreeHeta;

//                    var isBoundary = (i == 0 || i == NumberOfControlPointsKsi - 1 || j == 0 ||
//                                      j == NumberOfControlPointsHeta - 1);
//					collocationPoints.Add( new CollocationPoint2D(index++,coordinateKsi, coordinateHeta, isBoundary));
//				}
//			}
//#endregion

//            #region Knots
//            Vector singleKnotValuesKsi = KnotValueVectorKsi.RemoveDuplicatesFindMultiplicity()[0];
//			Vector singleKnotValuesHeta = KnotValueVectorHeta.RemoveDuplicatesFindMultiplicity()[0];

//			List<Knot> knots = new List<Knot>();

//			int id = 0;
//			for (int i = 0; i < singleKnotValuesKsi.Length; i++)
//			{
//				for (int j = 0; j < singleKnotValuesHeta.Length; j++)
//				{
//					knots.Add(new Knot() { ID = id, Ksi = singleKnotValuesKsi[i], Heta = singleKnotValuesHeta[j], Zeta = 0.0 });
//					id++;
//				}
//			}
//			#endregion

//			#region Elements
//			Vector multiplicityKsi = KnotValueVectorKsi.RemoveDuplicatesFindMultiplicity()[1];
//			Vector multiplicityHeta = KnotValueVectorHeta.RemoveDuplicatesFindMultiplicity()[1];

//			int numberOfElementsKsi = singleKnotValuesKsi.Length - 1;
//			int numberOfElementsHeta = singleKnotValuesHeta.Length - 1;
//			if (numberOfElementsKsi * numberOfElementsHeta == 0)
//			{
//				throw new NullReferenceException("Number of Elements should be defined before Element Connectivity");
//			}

//			for (int i = 0; i < numberOfElementsKsi; i++)
//			{
//				for (int j = 0; j < numberOfElementsHeta; j++)
//				{
//					IList<Knot> knotsOfElement = new List<Knot>();
//					knotsOfElement.Add(knots[i * singleKnotValuesHeta.Length + j]);
//					knotsOfElement.Add(knots[i * singleKnotValuesHeta.Length + j + 1]);
//					knotsOfElement.Add(knots[(i + 1) * singleKnotValuesHeta.Length + j]);
//					knotsOfElement.Add(knots[(i + 1) * singleKnotValuesHeta.Length + j + 1]);

//					int multiplicityElementKsi = 0;
//					if (multiplicityKsi[i + 1] - this.DegreeKsi > 0)
//					{
//						multiplicityElementKsi = (int)multiplicityKsi[i + 1] - this.DegreeKsi;
//					}

//					int multiplicityElementHeta = 0;
//					if (multiplicityHeta[j + 1] - this.DegreeHeta > 0)
//					{
//						multiplicityElementHeta = (int)multiplicityHeta[j + 1] - this.DegreeHeta;
//					}

//					int nurbsSupportKsi = this.DegreeKsi + 1;
//					int nurbsSupportHeta = this.DegreeHeta + 1;

//					IList<ControlPoint> elementControlPoints = new List<ControlPoint>();

//					for (int k = 0; k < nurbsSupportKsi; k++)
//					{
//						for (int l = 0; l < nurbsSupportHeta; l++)
//						{
//							int controlPointID = (i + multiplicityElementKsi) * this.NumberOfControlPointsHeta +
//								(j + multiplicityElementHeta) + k * this.NumberOfControlPointsHeta + l;
//							elementControlPoints.Add(this.ControlPoints[controlPointID]);
//						}
//					}

//					var elementCollocationPoints = collocationPoints.Where(cp =>
//						knotsOfElement[0].Ksi <= cp.Xi && cp.Xi <= knotsOfElement[3].Ksi &&
//						knotsOfElement[0].Heta <= cp.Eta && cp.Eta <= knotsOfElement[3].Heta).ToList();

//					foreach (var point in elementCollocationPoints)
//                    {
//                        if (Elements.Any(e => e.ID == point.ID)) continue;
//                        Element element = new NURBSElement2DCollocation
//                        {
//                            ID = point.ID,
//                            Patch = this,
//                            ElementType = new NURBSElement2DCollocation(),
//                            CollocationPoint = point
//                        };
//                        element.AddKnots(knotsOfElement);
//                        element.AddControlPoints(elementControlPoints.ToList<ControlPoint>());
//                        Elements.Add(element);
//                    }
//				}
//			}

//            var orderedElements=Elements.OrderBy(e => ((ICollocationElement) e).CollocationPoint.ID).ToList();
//            Elements.Clear();
//            Elements.AddRange(orderedElements);


//            #endregion
//        }

		private void CreatePatchData2D()
		{
			CreateNURBSElements2D();
			BuildEdgesDictionary();
		}

		private void CreateNURBSElements2D()
		{
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
					int elementID = i * numberOfElementsHeta + j;
					Element element = new NURBSElement2D()
					{
						ID = elementID,
						Patch = this,
						ElementType = new NURBSElement2D()
					};
					element.AddKnots(knotsOfElement);
					element.AddControlPoints(elementControlPoints.ToList<ControlPoint>());
					Elements.Add(element);
					//this.PatchesDictionary[1].ElementsDictionary.Add(element.ID, element);
				}
			}
			#endregion
		}

		private void CreateNURBSShells()
		{
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
							elementControlPoints.Add(ControlPoints[controlPointID]);
						}
					}
					int elementID = i * numberOfElementsHeta + j;
					Element element = new NURBSKirchhoffLoveShellElement()
					{
						ID = elementID,
						Patch = this,
						ElementType = new NURBSKirchhoffLoveShellElement()
					};
					element.AddKnots(knotsOfElement);
					element.AddControlPoints(elementControlPoints.ToList<ControlPoint>());
					Elements.Add(element);
					//this.PatchesDictionary[1].ElementsDictionary.Add(element.ID, element);
				}
			}
			#endregion
		}

		private void CreatePatchData3D()
		{
			CreateNURBSElements3D();
			BuildEdgesDictionary();
			BuildFacesDictionary();
		}

		private void CreateNURBSElements3D()
		{
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

						int elementID = i * numberOfElementsHeta * numberOfElementsZeta + j * numberOfElementsZeta + k;
						Element element = new NURBSElement3D()
						{
							ID = elementID,
							Patch = this,
							ElementType = new NURBSElement3D()
						};
						element.AddKnots(knotsOfElement);
						element.AddControlPoints(elementControlPoints.ToList<ControlPoint>());
						Elements.Add(element);
						//this.PatchesDictionary[0].ElementsDictionary.Add(element.ID, element);
					}
				}
			}
			#endregion
		}


		public void BuildFacesDictionary()
		{
			if (this.NumberOfDimensions <= 2)
			{
				FacesDictionary.Clear();
			}
			else
			{
				#region FaceRight
				Face faceRight = new Face();
				faceRight.Degrees[0] = this.DegreeHeta;
				faceRight.Degrees[1] = this.DegreeZeta;
				faceRight.KnotValueVectors.Add(0, KnotValueVectorHeta);
				faceRight.KnotValueVectors.Add(1, KnotValueVectorZeta);
				faceRight.Patch = this;
				int counter = 0;
				for (int i = 0; i < this.NumberOfControlPointsHeta; i++)
				{
					for (int j = 0; j < this.NumberOfControlPointsZeta; j++)
					{
						faceRight.ControlPointsDictionary.Add(counter++,
							ControlPoints[j + this.NumberOfControlPointsZeta * i]);
					}
				}
				FacesDictionary.Add(0, faceRight);
				#endregion
				#region FaceLeft
				Face faceLeft = new Face();
				faceLeft.Degrees[0] = this.DegreeHeta;
				faceLeft.Degrees[1] = this.DegreeZeta;
				faceLeft.KnotValueVectors.Add(0, KnotValueVectorHeta);
				faceLeft.KnotValueVectors.Add(1, KnotValueVectorZeta);
				faceLeft.Patch = this;
				counter = 0;
				for (int i = 0; i < this.NumberOfControlPointsHeta; i++)
				{
					for (int j = 0; j < this.NumberOfControlPointsZeta; j++)
					{
						faceLeft.ControlPointsDictionary.Add(counter++,
							ControlPoints[j + this.NumberOfControlPointsZeta * i +
							this.NumberOfControlPointsHeta * this.NumberOfControlPointsZeta * (this.NumberOfControlPointsKsi - 1)]);
					}
				}
				FacesDictionary.Add(1, faceLeft);
				#endregion
				#region FaceBottom
				Face faceBottom = new Face();
				faceBottom.Degrees[0] = this.DegreeKsi;
				faceBottom.Degrees[1] = this.DegreeHeta;
				faceBottom.KnotValueVectors.Add(0, KnotValueVectorKsi);
				faceBottom.KnotValueVectors.Add(1, KnotValueVectorHeta);
				faceBottom.Patch = this;
				counter = 0;
				for (int i = 0; i < this.NumberOfControlPointsKsi; i++)
				{
					for (int j = 0; j < this.NumberOfControlPointsHeta; j++)
					{
						faceBottom.ControlPointsDictionary.Add(counter++,
							ControlPoints[j * this.NumberOfControlPointsZeta +
							i * this.NumberOfControlPointsZeta * this.NumberOfControlPointsHeta]);
					}
				}
				FacesDictionary.Add(2, faceBottom);
				#endregion
				#region FaceUp
				Face faceUp = new Face();
				faceUp.Degrees[0] = this.DegreeKsi;
				faceUp.Degrees[1] = this.DegreeHeta;
				faceUp.KnotValueVectors.Add(0, KnotValueVectorKsi);
				faceUp.KnotValueVectors.Add(1, KnotValueVectorHeta);
				faceUp.Patch = this;
				counter = 0;
				for (int i = 0; i < this.NumberOfControlPointsKsi; i++)
				{
					for (int j = 0; j < this.NumberOfControlPointsHeta; j++)
					{
						faceUp.ControlPointsDictionary.Add(counter++,
							ControlPoints[this.NumberOfControlPointsZeta - 1 + j * this.NumberOfControlPointsZeta +
							i * this.NumberOfControlPointsZeta * this.NumberOfControlPointsHeta]);
					}
				}
				FacesDictionary.Add(3, faceUp);
				#endregion
				#region FaceFront
				Face faceFront = new Face();
				faceFront.Degrees[0] = this.DegreeKsi;
				faceFront.Degrees[1] = this.DegreeZeta;
				faceFront.KnotValueVectors.Add(0, KnotValueVectorKsi);
				faceFront.KnotValueVectors.Add(1, KnotValueVectorZeta);
				faceFront.Patch = this;
				counter = 0;
				for (int i = 0; i < this.NumberOfControlPointsKsi; i++)
				{
					for (int j = 0; j < this.NumberOfControlPointsZeta; j++)
					{
						faceFront.ControlPointsDictionary.Add(counter++,
							ControlPoints[j + i * this.NumberOfControlPointsHeta * this.NumberOfControlPointsZeta]);
					}
				}
				FacesDictionary.Add(4, faceFront);
				#endregion
				#region FaceBack
				Face faceBack = new Face();
				faceBack.Degrees[0] = this.DegreeKsi;
				faceBack.Degrees[1] = this.DegreeZeta;
				faceBack.KnotValueVectors.Add(0, KnotValueVectorKsi);
				faceBack.KnotValueVectors.Add(1, KnotValueVectorZeta);
				faceBack.Patch = this;
				counter = 0;
				for (int i = 0; i < this.NumberOfControlPointsKsi; i++)
				{
					for (int j = 0; j < this.NumberOfControlPointsZeta; j++)
					{
						faceBack.ControlPointsDictionary.Add(counter++,
							ControlPoints[j + i * this.NumberOfControlPointsHeta * this.NumberOfControlPointsZeta +
							this.NumberOfControlPointsZeta * (this.NumberOfControlPointsHeta - 1)]);
					}
				}
				FacesDictionary.Add(5, faceBack);
				#endregion
			}
		}

		internal void CreateNurbsShell()
		{
			BuildEdgesDictionary();
			CreateNURBSShells();
		}

		public void BuildEdgesDictionary()
		{
			if (this.NumberOfDimensions == 2)
			{
				#region EdgeRight
				Edge edgeRight = new Edge();
				edgeRight.ID = 0;
				edgeRight.Degree = this.DegreeHeta;
				edgeRight.KnotValueVector = this.KnotValueVectorHeta;
				edgeRight.numberOfControlPoints = this.NumberOfControlPointsHeta;
				edgeRight.Patch = this;
				int counter = 0;
				for (int i = 0; i < this.NumberOfControlPointsHeta; i++)
				{
					edgeRight.ControlPointsDictionary.Add(counter++, ControlPoints[i]);
				}
				EdgesDictionary.Add(0, edgeRight);
				#endregion

				#region EdgeLeft
				Edge edgeLeft = new Edge();
				edgeLeft.ID = 1;
				edgeLeft.Degree = this.DegreeHeta;
				edgeLeft.KnotValueVector = this.KnotValueVectorHeta;
				edgeLeft.numberOfControlPoints = this.NumberOfControlPointsHeta;
				edgeLeft.Patch = this;
				counter = 0;
				for (int i = 0; i < this.NumberOfControlPointsHeta; i++)
				{
					edgeLeft.ControlPointsDictionary.Add(counter++, ControlPoints[i + this.NumberOfControlPointsHeta * (this.NumberOfControlPointsKsi - 1)]);
				}
				EdgesDictionary.Add(1, edgeLeft);
				#endregion

				#region EdgeBottom
				Edge edgeBottom = new Edge();
				edgeBottom.ID = 2;
				edgeBottom.Degree = this.DegreeKsi;
				edgeBottom.KnotValueVector = this.KnotValueVectorKsi;
				edgeBottom.numberOfControlPoints = this.NumberOfControlPointsKsi;
				edgeBottom.Patch = this;
				counter = 0;
				for (int i = 0; i < this.NumberOfControlPointsKsi; i++)
				{
					edgeBottom.ControlPointsDictionary.Add(counter++, ControlPoints[i * this.NumberOfControlPointsHeta]);
				}
				EdgesDictionary.Add(2, edgeBottom);
				#endregion

				#region EdgeUp
				Edge edgeUp = new Edge();
				edgeUp.ID = 3;
				edgeUp.Degree = this.DegreeKsi;
				edgeUp.KnotValueVector = this.KnotValueVectorKsi;
				edgeUp.numberOfControlPoints = this.NumberOfControlPointsKsi;
				edgeUp.Patch = this;
				counter = 0;
				for (int i = 0; i < this.NumberOfControlPointsKsi; i++)
				{
					edgeUp.ControlPointsDictionary.Add(counter++,
						ControlPoints[i * this.NumberOfControlPointsHeta + this.NumberOfControlPointsHeta - 1]);
				}
				EdgesDictionary.Add(3, edgeUp);
				#endregion
			}
			else
			{
				#region Edge1
				Edge edge1 = new Edge();
				edge1.ID = 0;
				edge1.Degree = this.DegreeZeta;
				edge1.KnotValueVector = this.KnotValueVectorZeta;
				edge1.numberOfControlPoints = this.NumberOfControlPointsZeta;
				edge1.Patch = this;
				int counter = 0;
				for (int i = 0; i < this.NumberOfControlPointsZeta; i++)
				{
					edge1.ControlPointsDictionary.Add(counter++, ControlPoints[i]);
				}
				EdgesDictionary.Add(0, edge1);
				#endregion
				#region Edge2
				Edge edge2 = new Edge();
				edge2.ID = 1;
				edge2.Degree = this.DegreeZeta;
				edge2.KnotValueVector = this.KnotValueVectorZeta;
				edge2.numberOfControlPoints = this.NumberOfControlPointsZeta;
				edge2.Patch = this;
				counter = 0;
				for (int i = 0; i < this.NumberOfControlPointsZeta; i++)
				{
					edge2.ControlPointsDictionary.Add(counter++,
						ControlPoints[i + this.NumberOfControlPointsZeta * (this.NumberOfControlPointsHeta - 1)]);
				}
				EdgesDictionary.Add(1, edge2);
				#endregion
				#region Edge3
				Edge edge3 = new Edge();
				edge3.ID = 2;
				edge3.Degree = this.DegreeHeta;
				edge3.KnotValueVector = this.KnotValueVectorHeta;
				edge3.numberOfControlPoints = this.NumberOfControlPointsHeta;
				edge3.Patch = this;
				counter = 0;
				for (int i = 0; i < this.NumberOfControlPointsHeta; i++)
				{
					edge3.ControlPointsDictionary.Add(counter++,
						ControlPoints[i * this.NumberOfControlPointsZeta]);
				}
				EdgesDictionary.Add(2, edge3);
				#endregion
				#region Edge4
				Edge edge4 = new Edge();
				edge4.ID = 3;
				edge4.Degree = this.DegreeHeta;
				edge4.KnotValueVector = this.KnotValueVectorHeta;
				edge4.numberOfControlPoints = this.NumberOfControlPointsHeta;
				edge4.Patch = this;
				counter = 0;
				for (int i = 0; i < this.NumberOfControlPointsHeta; i++)
				{
					edge4.ControlPointsDictionary.Add(counter++,
						ControlPoints[i * this.NumberOfControlPointsZeta + this.NumberOfControlPointsZeta - 1]);
				}
				EdgesDictionary.Add(3, edge4);
				#endregion

				#region Edge5
				Edge edge5 = new Edge();
				edge5.ID = 4;
				edge5.Degree = this.DegreeZeta;
				edge5.KnotValueVector = this.KnotValueVectorZeta;
				edge5.numberOfControlPoints = this.NumberOfControlPointsZeta;
				int offset = this.NumberOfControlPointsZeta * this.NumberOfControlPointsHeta *
					(this.NumberOfControlPointsKsi - 1);
				edge5.Patch = this;
				counter = 0;
				for (int i = 0; i < this.NumberOfControlPointsZeta; i++)
				{
					edge5.ControlPointsDictionary.Add(counter++, ControlPoints[i + offset]);
				}
				EdgesDictionary.Add(4, edge5);
				#endregion
				#region Edge6
				Edge edge6 = new Edge();
				edge6.ID = 5;
				edge6.Degree = this.DegreeZeta;
				edge6.KnotValueVector = this.KnotValueVectorZeta;
				edge6.numberOfControlPoints = this.NumberOfControlPointsZeta;
				edge6.Patch = this;
				counter = 0;
				for (int i = 0; i < this.NumberOfControlPointsZeta; i++)
				{
					edge6.ControlPointsDictionary.Add(counter++,
						ControlPoints[i + this.NumberOfControlPointsZeta * (this.NumberOfControlPointsHeta - 1) + offset]);
				}
				EdgesDictionary.Add(5, edge6);
				#endregion
				#region Edge7
				Edge edge7 = new Edge();
				edge7.ID = 6;
				edge7.Degree = this.DegreeHeta;
				edge7.KnotValueVector = this.KnotValueVectorHeta;
				edge7.numberOfControlPoints = this.NumberOfControlPointsHeta;
				edge7.Patch = this;
				counter = 0;
				for (int i = 0; i < this.NumberOfControlPointsHeta; i++)
				{
					edge7.ControlPointsDictionary.Add(counter++,
						ControlPoints[i * this.NumberOfControlPointsZeta + offset]);
				}
				EdgesDictionary.Add(6, edge7);
				#endregion
				#region Edge8
				Edge edge8 = new Edge();
				edge8.ID = 7;
				edge8.Degree = this.DegreeHeta;
				edge8.KnotValueVector = this.KnotValueVectorHeta;
				edge8.numberOfControlPoints = this.NumberOfControlPointsHeta;
				edge8.Patch = this;
				counter = 0;
				for (int i = 0; i < this.NumberOfControlPointsHeta; i++)
				{
					edge8.ControlPointsDictionary.Add(counter++,
						ControlPoints[i * this.NumberOfControlPointsZeta + this.NumberOfControlPointsZeta - 1 + offset]);
				}
				EdgesDictionary.Add(7, edge8);
				#endregion

				#region Edge9
				Edge edge9 = new Edge();
				edge9.ID = 8;
				edge9.Degree = this.DegreeKsi;
				edge9.KnotValueVector = this.KnotValueVectorKsi;
				edge9.numberOfControlPoints = this.NumberOfControlPointsKsi;
				edge9.Patch = this;
				counter = 0;
				for (int i = 0; i < this.NumberOfControlPointsKsi; i++)
				{
					edge9.ControlPointsDictionary.Add(counter++,
						ControlPoints[i * this.NumberOfControlPointsZeta * this.NumberOfControlPointsHeta]);
				}
				EdgesDictionary.Add(8, edge9);
				#endregion
				#region Edge10
				Edge edge10 = new Edge();
				edge10.ID = 9;
				edge10.Degree = this.DegreeKsi;
				edge10.KnotValueVector = this.KnotValueVectorKsi;
				edge10.numberOfControlPoints = this.NumberOfControlPointsKsi;
				edge10.Patch = this;
				counter = 0;
				for (int i = 0; i < this.NumberOfControlPointsKsi; i++)
				{
					edge10.ControlPointsDictionary.Add(counter++,
						ControlPoints[i * this.NumberOfControlPointsZeta * this.NumberOfControlPointsHeta
						+ this.NumberOfControlPointsZeta - 1]);
				}
				EdgesDictionary.Add(9, edge10);
				#endregion
				#region Edge11
				Edge edge11 = new Edge();
				edge11.ID = 10;
				edge11.Degree = this.DegreeKsi;
				edge11.KnotValueVector = this.KnotValueVectorKsi;
				edge11.numberOfControlPoints = this.NumberOfControlPointsKsi;
				edge11.Patch = this;
				counter = 0;
				for (int i = 0; i < this.NumberOfControlPointsKsi; i++)
				{
					edge11.ControlPointsDictionary.Add(counter++,
						ControlPoints[i * this.NumberOfControlPointsZeta * this.NumberOfControlPointsHeta +
						this.NumberOfControlPointsZeta * (this.NumberOfControlPointsHeta - 1)]);
				}
				EdgesDictionary.Add(10, edge11);
				#endregion
				#region Edge12
				Edge edge12 = new Edge();
				edge12.ID = 11;
				edge12.Degree = this.DegreeKsi;
				edge12.KnotValueVector = this.KnotValueVectorKsi;
				edge12.numberOfControlPoints = this.NumberOfControlPointsKsi;
				edge12.Patch = this;
				counter = 0;
				for (int i = 0; i < this.NumberOfControlPointsKsi; i++)
				{
					edge12.ControlPointsDictionary.Add(counter++,
						ControlPoints[i * this.NumberOfControlPointsZeta * this.NumberOfControlPointsHeta +
						this.NumberOfControlPointsZeta * (this.NumberOfControlPointsHeta - 1) + this.NumberOfControlPointsZeta - 1]);
				}
				EdgesDictionary.Add(11, edge12);
				#endregion
			}

		}
	}
}
