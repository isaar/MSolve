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

namespace ISAAR.MSolve.IGA.Entities
{
    public class Face : Boundary
    {
        public const int NumberOfDimensions =2;

        public int numberOfControlPoints { get; set; }

        public int[] Degrees = new int[2];

        public Patch Patch { get; set; }

        public Dictionary<int,IVector> KnotValueVectors =new Dictionary<int, IVector>();

        private readonly Dictionary<int, ControlPoint> controlPointsDictionary = new Dictionary<int, ControlPoint>();

        private readonly Dictionary<int, Element> elementsDictionary = new Dictionary<int, Element>();

        private readonly Dictionary< int , Dictionary<DOFType,int>> controlPointsDOFsDictionary = new Dictionary<int, Dictionary<DOFType, int>>();

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
            Dictionary<int, double> faceLoad = new Dictionary<int, double>();
            foreach (LoadingCondition loading in loadingConditions)
            {
                Dictionary<int, double> load = CalculateLoadingCondition(loading);
                foreach (int dof in load.Keys)
                {
                    if (faceLoad.ContainsKey(dof))
                    {
                        faceLoad[dof] += load[dof];
                    }
                    else
                    {
                        faceLoad.Add(dof, load[dof]);
                    }
                }
            }
            return faceLoad;
        }

        private Dictionary<int, double> CalculateLoadingCondition(LoadingCondition loading)
        {
            LoadProvider provider = new LoadProvider();
            if (elementsDictionary.Count == 0) CreateFaceElements();
            Dictionary<int, double> load = new Dictionary<int, double>();
            if (loading is NeumannBoundaryCondition)
            {
                foreach (Element element in elementsDictionary.Values)
                    foreach (int dof in provider.LoadNeumann(element, this, loading as NeumannBoundaryCondition).Keys)
                    {
                        if (load.ContainsKey(dof))
                        {
                            load[dof] += provider.LoadNeumann(element, this, loading as NeumannBoundaryCondition)[dof];
                        }
                        else
                        {
                            load.Add(dof, provider.LoadNeumann(element, this, loading as NeumannBoundaryCondition)[dof]);
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

	    private void CreateFaceElements()
        {
            #region Knots
            Vector singleKnotValuesKsi = KnotValueVectors[0].RemoveDuplicatesFindMultiplicity()[0];
            Vector singleKnotValuesHeta = KnotValueVectors[1].RemoveDuplicatesFindMultiplicity()[0];

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
            Vector multiplicityKsi = KnotValueVectors[0].RemoveDuplicatesFindMultiplicity()[1];
            Vector multiplicityHeta = KnotValueVectors[1].RemoveDuplicatesFindMultiplicity()[1];

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
                    if (multiplicityKsi[i + 1] - this.Degrees[0] > 0)
                    {
                        multiplicityElementKsi = (int)multiplicityKsi[i + 1] - this.Degrees[0];
                    }

                    int multiplicityElementHeta = 0;
                    if (multiplicityHeta[j + 1] - this.Degrees[1] > 0)
                    {
                        multiplicityElementHeta = (int)multiplicityHeta[j + 1] - this.Degrees[1];
                    }

                    int nurbsSupportKsi = this.Degrees[0] + 1;
                    int nurbsSupportHeta = this.Degrees[1] + 1;

                    int NumberOfControlPointsHeta = KnotValueVectors[1].Length - Degrees[1] - 1;

                    IList<ControlPoint> elementControlPoints = new List<ControlPoint>();

                    for (int k = 0; k < nurbsSupportKsi; k++)
                    {
                        for (int l = 0; l < nurbsSupportHeta; l++)
                        {
                            int controlPointID = (i + multiplicityElementKsi) * NumberOfControlPointsHeta +
                                (j + multiplicityElementHeta) + k * NumberOfControlPointsHeta + l;

							elementControlPoints.Add(this.controlPointsDictionary[controlPointID]);
                        }
                    }
                    int elementID = i * numberOfElementsHeta + j;
                    Element element = new NURBSElement2D()
                    {
                        ID = elementID,
                        ElementType = new NURBSElement2D(),
						Patch = Patch
                    };
                    element.AddKnots(knotsOfElement);
                    element.AddControlPoints(elementControlPoints);
                    this.elementsDictionary.Add(elementID, element);
                    //this.PatchesDictionary[1].ElementsDictionary.Add(element.ID, element);
                }
            }
            #endregion
        }

    }
}
