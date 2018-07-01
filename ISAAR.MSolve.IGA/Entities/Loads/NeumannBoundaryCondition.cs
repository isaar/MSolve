using ISAAR.MSolve.IGA.Problems.SupportiveClasses;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.IGA.Entities.Loads
{
    public class NeumannBoundaryCondition : LoadingCondition
    {
        public Value Value { get; private set; }
        
        public NeumannBoundaryCondition(Value neumannValue)
        {
            this.Value = neumannValue;
        }

    }

    //public delegate double Value(GaussLegendrePoint3D gaussPoint);

    public delegate double[] Value(double x, double y, double z);

    public class LoadProvider
    {
        public Dictionary<int, double> LoadNeumann(Element element,Edge edge, NeumannBoundaryCondition neumann)
        {
            return element.ElementType.CalculateLoadingCondition(element,edge, neumann);
        }

        public Dictionary<int, double> LoadNeumann(Element element, Face face, NeumannBoundaryCondition neumann)
        {
            return element.ElementType.CalculateLoadingCondition(element, face, neumann);
        }

        public Dictionary<int, double> LoadPressure(Element element, Edge edge, PressureBoundaryCondition pressure)
        {
            return element.ElementType.CalculateLoadingCondition(element, edge, pressure);
        }

        public Dictionary<int, double> LoadPressure(Element element, Face face, PressureBoundaryCondition pressure)
        {
            return element.ElementType.CalculateLoadingCondition(element, face, pressure);
        }
    }
}
