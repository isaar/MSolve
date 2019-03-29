using ISAAR.MSolve.IGA.Entities.Loads;
using ISAAR.MSolve.IGA.Interfaces;
using ISAAR.MSolve.IGA.Problems.Structural.Elements;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.IGA.Elements;

namespace ISAAR.MSolve.IGA.Entities
{
	public class Edge : Boundary
	{
		public int ID { get; set; }

		public const int NumberOfDimensions = 1;

		public int numberOfControlPoints { get; set; }

		public int Degree { get; set; }

		public Patch Patch { get; set; }
		//public Patch_v2 Patch_v2 { get; set; }

		public IVector KnotValueVector { get; set; }

		private readonly Dictionary<int, ControlPoint> controlPointsDictionary = new Dictionary<int, ControlPoint>();

		private readonly Dictionary<int, Element> elementsDictionary = new Dictionary<int, Element>();

		private readonly Dictionary<int, Dictionary<DOFType, int>> controlPointsDOFsDictionary =
			new Dictionary<int, Dictionary<DOFType, int>>();

		private readonly List<IBoundaryCondition> boundaryConditions = new List<IBoundaryCondition>();

		private readonly List<LoadingCondition> loadingConditions = new List<LoadingCondition>();


		#region Properties

		public Dictionary<int, ControlPoint> ControlPointsDictionary
		{
			get { return controlPointsDictionary; }
		}

		public Dictionary<int, Element> ElementsDictionary
		{
			get { return elementsDictionary; }
		}

		public Dictionary<int, Dictionary<DOFType, int>> ControlPointDOFsDictionary
		{
			get { return controlPointsDOFsDictionary; }
		}

		public List<IBoundaryCondition> BoundaryConditions
		{
			get { return boundaryConditions; }
		}

		public List<LoadingCondition> LoadingConditions
		{
			get { return loadingConditions; }
		}

		#endregion

		public Dictionary<int, double> CalculateLoads()
		{
			Dictionary<int, double> edgeLoad = new Dictionary<int, double>();
			foreach (LoadingCondition loading in loadingConditions)
			{
				Dictionary<int, double> load = CalculateLoadingCondition(loading);
				foreach (int dof in load.Keys)
				{
					if (edgeLoad.ContainsKey(dof))
					{
						edgeLoad[dof] += load[dof];
					}
					else
					{
						edgeLoad.Add(dof, load[dof]);
					}
				}
			}

			return edgeLoad;
		}

		private Dictionary<int, double> CalculateLoadingCondition(LoadingCondition loading)
		{
			LoadProvider provider = new LoadProvider();
			if (elementsDictionary.Count == 0) CreateEdgeElements();
			Dictionary<int, double> load = new Dictionary<int, double>();
			if (loading is NeumannBoundaryCondition)
			{
				foreach (Element element in elementsDictionary.Values)
				{
					var loadNeumann = provider.LoadNeumann(element, this, loading as NeumannBoundaryCondition);
					foreach (int dof in loadNeumann.Keys)
					{
						if (load.ContainsKey(dof))
						{
							load[dof] += loadNeumann[dof];
						}
						else
						{
							load.Add(dof, loadNeumann[dof]);
						}
					}
				}
			}
			else if (loading is PressureBoundaryCondition)
			{
				foreach (Element element in elementsDictionary.Values)
				foreach (int dof in provider.LoadPressure(element, this, loading as PressureBoundaryCondition).Keys)
				{
					if (load.ContainsKey(dof))
					{
						load[dof] += provider.LoadPressure(element, this, loading as PressureBoundaryCondition)[dof];
					}
					else
					{
						load.Add(dof, provider.LoadPressure(element, this, loading as PressureBoundaryCondition)[dof]);
					}
				}
			}

			return load;
		}

		private void CreateEdgeElements()
		{
			#region Knots

			Vector singleKnotValuesKsi = KnotValueVector.RemoveDuplicatesFindMultiplicity()[0];

			List<Knot> knots = new List<Knot>();

			int id = 0;
			for (int i = 0; i < singleKnotValuesKsi.Length; i++)
			{
				knots.Add(new Knot() {ID = id, Ksi = singleKnotValuesKsi[i], Heta = 0.0, Zeta = 0.0});
				id++;
			}

			#endregion

			#region Elements

			Vector multiplicityKsi = KnotValueVector.RemoveDuplicatesFindMultiplicity()[1];

			int numberOfElementsKsi = singleKnotValuesKsi.Length - 1;
			if (numberOfElementsKsi == 0)
			{
				throw new NullReferenceException("Number of Elements should be defined before Element Connectivity");
			}

			for (int i = 0; i < numberOfElementsKsi; i++)
			{
				IList<Knot> knotsOfElement = new List<Knot>();
				knotsOfElement.Add(knots[i]);
				knotsOfElement.Add(knots[i + 1]);

				int multiplicityElementKsi = 0;
				if (multiplicityKsi[i + 1] - this.Degree > 0)
				{
					multiplicityElementKsi = (int) multiplicityKsi[i + 1] - this.Degree;
				}

				int nurbsSupportKsi = this.Degree + 1;

				IList<ControlPoint> elementControlPoints = new List<ControlPoint>();

				for (int k = 0; k < nurbsSupportKsi; k++)
				{
					int controlPointID = i + multiplicityElementKsi + k;
					elementControlPoints.Add(this.controlPointsDictionary[controlPointID]);
				}

				int elementID = i;
				Element element = new NURBSElement1D()
				{
					ID = elementID,
					ElementType = new NURBSElement1D(),
					Patch = this.Patch,
					Degree = this.Degree,
					Model = Patch.Elements[0].Model
				};

				element.AddKnots(knotsOfElement);
				element.AddControlPoints(elementControlPoints);
				this.elementsDictionary.Add(elementID, element);
				//this.PatchesDictionary[1].ElementsDictionary.Add(element.ID, element);
			}

			#endregion
		}
	}
}