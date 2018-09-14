using ISAAR.MSolve.IGA.Entities;
using ISAAR.MSolve.IGA.Interfaces;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.IGA.Entities.Loads;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using ISAAR.MSolve.IGA.Problems.SupportiveClasses;
using ISAAR.MSolve.Numerical.LinearAlgebra;

namespace ISAAR.MSolve.IGA.Elements
{
    public class NURBSElement1D : Element, IStructuralIsogeometricElement
    {

        protected readonly static DOFType[] controlPointDOFTypes = new DOFType[] { DOFType.X};
        protected DOFType[][] dofTypes;
        protected IElementDOFEnumerator dofEnumerator = new GenericDOFEnumerator();
        public  int Degree { get; set; }

        #region IStructuralIsogeometricElement
        public ElementDimensions ElementDimensions { get { return ElementDimensions.OneD; } }

        public IElementDOFEnumerator DOFEnumerator { get { return dofEnumerator; } set { this.dofEnumerator = value; } }

        public bool MaterialModified => throw new NotImplementedException();

        public double[] CalculateForces(Element element, double[] localDisplacements, double[] localdDisplacements)
        {
            throw new NotImplementedException();
        }

        public double[] CalculateForcesForLogging(Element element, double[] localDisplacements)
        {
            throw new NotImplementedException();
        }

        public Dictionary<int, double> CalculateLoadingCondition(Element element,Edge edge, NeumannBoundaryCondition neumann)
        {
            IList<GaussLegendrePoint3D> gaussPoints = CreateElementGaussPoints(element);
            Dictionary<int, double> neumannLoad = new Dictionary<int, double>();
            IList<ControlPoint> controlPoints = new List<ControlPoint>();

            foreach (ControlPoint controlPoint in element.ControlPoints)
            {
                if (element.Patch.NumberOfDimensions == 2)
                {
                    controlPoints.Add(new ControlPoint()
                    {
                        ID = (edge.ID < 2) ? controlPoint.ID % element.Patch.NumberOfControlPointsHeta: controlPoint.ID / element.Patch.NumberOfControlPointsHeta ,
                        Ksi = controlPoint.Ksi,
                        Heta = controlPoint.Heta,
                        Zeta = controlPoint.Zeta,
                        X = controlPoint.X,
                        Y = controlPoint.Y,
                        Z = controlPoint.Z,
                        WeightFactor = controlPoint.WeightFactor                        
                    });
                }
                else
                {
                    int ID=-1;
                    switch (edge.ID)
                    {
                        case 1:
                        case 2:
                        case 5:
                        case 6:
                            ID= controlPoint.ID % (element.Patch.NumberOfControlPointsHeta * element.Patch.NumberOfControlPointsZeta) % element.Patch.NumberOfControlPointsZeta;
                            break;
                        case 3:
                        case 4:
                        case 7:
                        case 8:
                            ID= controlPoint.ID % (element.Patch.NumberOfControlPointsHeta * element.Patch.NumberOfControlPointsZeta) / element.Patch.NumberOfControlPointsZeta;
                            break;

                        case 9:
                        case 10:
                        case 11:
                        case 12:
                            ID = controlPoint.ID / (element.Patch.NumberOfControlPointsHeta * element.Patch.NumberOfControlPointsZeta);
                            break;
                    }
                    controlPoints.Add(new ControlPoint()
                    {
                        ID = ID,
                        Ksi = controlPoint.Ksi,
                        Heta = controlPoint.Heta,
                        Zeta = controlPoint.Zeta,
                        X = controlPoint.X,
                        Y = controlPoint.Y,
                        Z = controlPoint.Z,
                        WeightFactor = controlPoint.WeightFactor
                    });

                }
            }


            NURBS1D nurbs = new NURBS1D(element, controlPoints, edge);

            for (int j = 0; j < gaussPoints.Count; j++)
            {
                double xGaussPoint = 0;
                double yGaussPoint = 0;
                double zGaussPoint = 0;
                double jacobian1 = 0.0;
                double jacobian2 = 0.0;
                for (int k = 0; k < element.ControlPoints.Count; k++)
                {
                    xGaussPoint += nurbs.NurbsValues[k, j] * element.ControlPoints[k].X;
                    yGaussPoint += nurbs.NurbsValues[k, j] * element.ControlPoints[k].Y;
                    zGaussPoint += nurbs.NurbsValues[k, j] * element.ControlPoints[k].Z;
                    jacobian1 += nurbs.NurbsDerivativeValuesKsi[k, j] * element.ControlPoints[k].X;
                    jacobian2 += nurbs.NurbsDerivativeValuesKsi[k, j] * element.ControlPoints[k].Y;
                }
                double jacdet = Math.Sqrt(Math.Pow(jacobian1, 2) + Math.Pow(jacobian2, 2));
                var loadGaussPoint = neumann.Value(xGaussPoint, yGaussPoint, zGaussPoint);

                for (int k = 0; k < element.ControlPoints.Count; k++)
                {
                    int dofIDX = element.Patch.ControlPointDOFsDictionary[element.ControlPoints[k].ID][DOFType.X];
                    int dofIDY = element.Patch.ControlPointDOFsDictionary[element.ControlPoints[k].ID][DOFType.Y];
                    if (neumannLoad.ContainsKey(dofIDX))
                        neumannLoad[dofIDX] += jacdet * gaussPoints[j].WeightFactor * nurbs.NurbsValues[k, j] * loadGaussPoint[0];
                    else
                        neumannLoad.Add(dofIDX, jacdet * gaussPoints[j].WeightFactor * nurbs.NurbsValues[k, j] * loadGaussPoint[0]);

                    if (neumannLoad.ContainsKey(dofIDY))
                        neumannLoad[dofIDY] += jacdet * gaussPoints[j].WeightFactor * nurbs.NurbsValues[k, j] * loadGaussPoint[1];
                    else
                        neumannLoad.Add(dofIDY, jacdet * gaussPoints[j].WeightFactor * nurbs.NurbsValues[k, j] * loadGaussPoint[1]);
                }


            }
            return neumannLoad;
        }


        public Dictionary<int, double> CalculateLoadingCondition(Element element, Edge edge, PressureBoundaryCondition pressure)
        {
            IList<GaussLegendrePoint3D> gaussPoints = CreateElementGaussPoints(element);
            Dictionary<int, double> pressureLoad = new Dictionary<int, double>();
            IList<ControlPoint> controlPoints = new List<ControlPoint>();

            foreach (ControlPoint controlPoint in element.ControlPoints)
            {
                if (element.Patch.NumberOfDimensions == 2)
                {
                    controlPoints.Add(new ControlPoint()
                    {
                        ID = (edge.ID < 2) ? controlPoint.ID % element.Patch.NumberOfControlPointsHeta : controlPoint.ID / element.Patch.NumberOfControlPointsHeta,
                        Ksi = controlPoint.Ksi,
                        Heta = controlPoint.Heta,
                        Zeta = controlPoint.Zeta,
                        X = controlPoint.X,
                        Y = controlPoint.Y,
                        Z = controlPoint.Z,
                        WeightFactor = controlPoint.WeightFactor
                    });
                }
                else
                {
                    int ID = -1;
                    switch (edge.ID)
                    {
                        case 1:
                        case 2:
                        case 5:
                        case 6:
                            ID = controlPoint.ID % (element.Patch.NumberOfControlPointsHeta * element.Patch.NumberOfControlPointsZeta) % element.Patch.NumberOfControlPointsZeta;
                            break;
                        case 3:
                        case 4:
                        case 7:
                        case 8:
                            ID = controlPoint.ID % (element.Patch.NumberOfControlPointsHeta * element.Patch.NumberOfControlPointsZeta) / element.Patch.NumberOfControlPointsZeta;
                            break;

                        case 9:
                        case 10:
                        case 11:
                        case 12:
                            ID = controlPoint.ID / (element.Patch.NumberOfControlPointsHeta * element.Patch.NumberOfControlPointsZeta);
                            break;
                    }
                    controlPoints.Add(new ControlPoint()
                    {
                        ID = ID,
                        Ksi = controlPoint.Ksi,
                        Heta = controlPoint.Heta,
                        Zeta = controlPoint.Zeta,
                        X = controlPoint.X,
                        Y = controlPoint.Y,
                        Z = controlPoint.Z,
                        WeightFactor = controlPoint.WeightFactor
                    });

                }
            }


            NURBS1D nurbs = new NURBS1D(element, controlPoints, edge);

            for (int j = 0; j < gaussPoints.Count; j++)
            {
                double xGaussPoint = 0;
	            double yGaussPoint = 0;
                double jacobian1 = 0.0;
                double jacobian2 = 0.0;
                for (int k = 0; k < element.ControlPoints.Count; k++)
                {
                    xGaussPoint += nurbs.NurbsValues[k, j] * element.ControlPoints[k].X;
                    yGaussPoint += nurbs.NurbsValues[k, j] * element.ControlPoints[k].Y;
                    jacobian1 += nurbs.NurbsDerivativeValuesKsi[k, j] * element.ControlPoints[k].X;
                    jacobian2 += nurbs.NurbsDerivativeValuesKsi[k, j] * element.ControlPoints[k].Y;
                }
                double jacdet = Math.Sqrt(Math.Pow(jacobian1, 2) + Math.Pow(jacobian2, 2));

                double norm = Math.Sqrt(Math.Pow(xGaussPoint, 2) + Math.Pow(yGaussPoint, 2));
                var loadGaussPointX = pressure.Value* xGaussPoint / norm;
	            var loadGaussPointY = pressure.Value* yGaussPoint / norm;

				for (int k = 0; k < element.ControlPoints.Count; k++)
                {
                    int dofIDX = element.Patch.ControlPointDOFsDictionary[element.ControlPoints[k].ID][DOFType.X];
                    int dofIDY = element.Patch.ControlPointDOFsDictionary[element.ControlPoints[k].ID][DOFType.Y];
                    if (pressureLoad.ContainsKey(dofIDX))
                        pressureLoad[dofIDX] += jacdet * gaussPoints[j].WeightFactor * nurbs.NurbsValues[k, j] * loadGaussPointX;
                    else
                        pressureLoad.Add(dofIDX, jacdet * gaussPoints[j].WeightFactor * nurbs.NurbsValues[k, j] *  loadGaussPointX);

                    if (pressureLoad.ContainsKey(dofIDY))
                        pressureLoad[dofIDY] += jacdet * gaussPoints[j].WeightFactor * nurbs.NurbsValues[k, j] *  loadGaussPointY;
                    else
                        pressureLoad.Add(dofIDY, jacdet * gaussPoints[j].WeightFactor * nurbs.NurbsValues[k, j] *  loadGaussPointY);
                }
            }
            return pressureLoad;
        }
        
        private IList<GaussLegendrePoint3D> CreateElementGaussPoints(Element element)
        {
            GaussQuadrature gauss = new GaussQuadrature();
            return gauss.CalculateElementGaussPoints(((NURBSElement1D)element).Degree, element.Knots);
        }

        

        public Tuple<double[], double[]> CalculateStresses(Element element, double[] localDisplacements, double[] localdDisplacements)
        {
            throw new NotImplementedException();
        }

        public void ClearMaterialState()
        {
            throw new NotImplementedException();
        }

        public IMatrix2D DampingMatrix(Element element)
        {
            throw new NotImplementedException();
        }

        public IList<IList<DOFType>> GetElementDOFTypes(Element element)
        {
            dofTypes = new DOFType[element.ControlPoints.Count][];
            for (int i = 0; i < element.ControlPoints.Count; i++)
            {
                dofTypes[i] = controlPointDOFTypes;
            }
            return dofTypes;
        }

        public IMatrix2D MassMatrix(Element element)
        {
            throw new NotImplementedException();
        }

        public void ResetMaterialModified()
        {
            throw new NotImplementedException();
        }

        public IMatrix2D StiffnessMatrix(IElement element)
        {
            throw new NotImplementedException();
        }

        public Dictionary<int, double> CalculateLoadingCondition(Element element, Face face, NeumannBoundaryCondition neumann)
        {
            throw new NotSupportedException();
        }

        public Dictionary<int, double> CalculateLoadingCondition(Element element, Face face, PressureBoundaryCondition pressure)
        {
            throw new NotSupportedException();
        }

		public IMatrix2D StiffnessMatrix(Element element)
		{
			throw new NotImplementedException();
		}

		public IMatrix2D MassMatrix(IElement element)
		{
			throw new NotImplementedException();
		}

		public IMatrix2D DampingMatrix(IElement element)
		{
			throw new NotImplementedException();
		}

		public IList<IList<DOFType>> GetElementDOFTypes(IElement element)
		{
			throw new NotImplementedException();
		}
		#endregion
	}
}
