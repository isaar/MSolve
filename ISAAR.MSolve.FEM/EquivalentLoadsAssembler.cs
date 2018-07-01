using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.Numerical.LinearAlgebra;
using ISAAR.MSolve.Numerical.LinearAlgebra.Interfaces;
using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.FEM
{
    public class EquivalentLoadsAssembler
    {
        public static void AssignEquivalentNodalLoads(Subdomain subdomain, IElementMatrixProvider elementProvider, IVector solution, IVector dSolution) //TODOMaria this should also take as argument the nodal displacements of the constraints (after the refactoring)
        {
            var times = new Dictionary<string, TimeSpan>();
            var totalStart = DateTime.Now;
            times.Add("rowIndexCalculation", DateTime.Now - totalStart);
            times.Add("element", TimeSpan.Zero);
            times.Add("addition", TimeSpan.Zero);
            var subdomainForces = new double[subdomain.Forces.Length];
            foreach (Element element in subdomain.ElementsDictionary.Values)
            {
                var isEmbeddedElement = element.ElementType is IEmbeddedElement;
                var elStart = DateTime.Now;
                IMatrix2D ElementK = elementProvider.Matrix(element);

                double[] localSolution = subdomain.CalculateElementNodalDisplacements(element, solution);
                double[] localdSolution = subdomain.CalculateElementNodalDisplacements(element, dSolution);

                var equivalentNodalForces = new double[localSolution.Length];
                ElementK.Multiply(new Vector(localSolution), equivalentNodalForces);
                subdomain.AddLocalVectorToGlobal(element, equivalentNodalForces, subdomainForces);

                times["addition"] += DateTime.Now - elStart;
            }

            for (int i = 0; i < subdomainForces.Length; i++)//TODOMaria is this right?? maybe we should keep a separate copy of subdomain forces in each subdomain
                subdomain.Forces[i] += subdomainForces[i];

            var totalTime = DateTime.Now - totalStart;
        }
    }
}
