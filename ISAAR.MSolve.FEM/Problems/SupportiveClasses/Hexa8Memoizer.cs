using System;
using System.Collections.Generic;

namespace ISAAR.MSolve.FEM.Problems.SupportiveClasses
{
    public class Hexa8Memoizer
    {
        private readonly Dictionary<int, Tuple<double[], double[,,]>> integrationDictionary = new Dictionary<int, Tuple<double[], double[,,]>>();

        public Tuple<double[], double[,,]> GetIntegrationData(int element)
        {
            if (integrationDictionary.ContainsKey(element))
                return integrationDictionary[element];
            else
                return new Tuple<double[], double[,,]>(null, null);
        }

        public void SetIntegrationData(int element, Tuple<double[], double[,,]> integrationData)
        {
            integrationDictionary.Add(element, integrationData);
        }
    }

}
