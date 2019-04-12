using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.LinearAlgebra.Output;

namespace ISAAR.MSolve.MultiscaleAnalysisMerge.SupportiveClasses
{
    /// <summary>
    /// Model separation methods for models with embedded elements
    /// Authors: Gerasimos Sotiropoulos
    /// </summary>
    public static class DdmCalculations // exei allaxei mono to onoma den exei ginei update se model xrhsh klp.
    {
        public static int[][] CalculateSubdElementIds(int hexa1, int hexa2, int hexa3, int elem1, int elem2, Model model)
        {
            int[][] subdElementIds = new int[8][];
            for (int i1 = 0; i1 < 7; i1++)
            {subdElementIds[i1] = new int[hexa1 * hexa2 * hexa3 / 8];}
            subdElementIds[7] = new int[(hexa1 * hexa2 * hexa3 / 8) + 3 * elem1 * elem2];

            int [] subdElementCounters = new int[8];

            for (int h1 = 0; h1 < hexa1; h1++)
            {
                for (int h2 = 0; h2 < hexa2; h2++)
                {
                    for (int h3 = 0; h3 < hexa3; h3++)
                    {
                        int ElementID = h1 + 1 + (h2 + 1 - 1) * hexa1 + (h3 + 1 - 1) * (hexa1) * hexa2; // h1+1 dioti h1 einai zero based

                        int s1; int s2; int s3;
                        if (h1 <= 0.5 * hexa1-1) { s1 = 1; } else { s1 = 2; };
                        if (h2 <= 0.5 * hexa2-1) { s2 = 1; } else { s2 = 2; };
                        if (h3 <= 0.5 * hexa3-1) { s3 = 1; } else { s3 = 2; };

                        int subdID = s1 + (s2 - 1) * 2 + (s3 - 1) * 4;

                        subdElementIds[subdID - 1][subdElementCounters[subdID - 1]] = ElementID;
                        subdElementCounters[subdID - 1] += 1;
                        //model.ElementsDictionary.Add(e1.ID, e1);
                        //model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(e1.ID, e1);

                    }
                }
            }

            
            //int subdID = 8;

            for (int ElementID = hexa1 * hexa2 * hexa3 + 1; ElementID < hexa1 * hexa2 * hexa3+ 3*(elem1*elem2)+1; ElementID++)
            {
                subdElementIds[7][subdElementCounters[7]] = ElementID;
                subdElementCounters[7] += 1;
            }
            return subdElementIds;
        }

        public static int[][] CalculateSubdElementIdsVerticalHexaOnly(int hexa1, int hexa2, int hexa3, int elem1, int elem2, Model model,int nSubdomains)
        {
            //int[][] subdElementIds = new int[8][];
            int[][] subdElementIds = new int[nSubdomains][];
            for (int i1 = 0; i1 < nSubdomains; i1++)
            { subdElementIds[i1] = new int[hexa1 * hexa2 * hexa3 / nSubdomains]; }
           


            int[] subdElementCounters = new int[8];

            for (int h1 = 0; h1 < hexa1; h1++)
            {
                for (int h2 = 0; h2 < hexa2; h2++)
                {
                    for (int h3 = 0; h3 < hexa3; h3++)
                    {
                        int ElementID = h1 + 1 + (h2 + 1 - 1) * hexa1 + (h3 + 1 - 1) * (hexa1) * hexa2; // h1+1 dioti h1 einai zero based

                        int s1; int s2; int s3;
                        if (h1 <= 0.5 * hexa1 - 1) { s1 = 1; } else { s1 = 2; };
                        if (h2 <= 0.5 * hexa2 - 1) { s2 = 1; } else { s2 = 2; };
                        if (h3 <= 0.5 * hexa3 - 1) { s3 = 1; } else { s3 = 2; };

                        int subdID = s1 + (s2 - 1) * 2 + (s3 - 1) * 4;

                        subdElementIds[subdID - 1][subdElementCounters[subdID - 1]] = ElementID;
                        subdElementCounters[subdID - 1] += 1;
                        //model.ElementsDictionary.Add(e1.ID, e1);
                        //model.SubdomainsDictionary[subdomainID].ElementsDictionary.Add(e1.ID, e1);

                    }
                }
            }


            //int subdID = 8;

            for (int ElementID = hexa1 * hexa2 * hexa3 + 1; ElementID < hexa1 * hexa2 * hexa3 + 3 * (elem1 * elem2) + 1; ElementID++)
            {
                subdElementIds[7][subdElementCounters[7]] = ElementID;
                subdElementCounters[7] += 1;
            }
            return subdElementIds;
        }

        public static void SeparateSubdomains(Model model, int[][] subdElementIds)
        {
            model.SubdomainsDictionary.Clear();

            for (int subdID = 0; subdID < subdElementIds.GetLength(0); subdID++)
            {
                model.SubdomainsDictionary.Add(subdID, new Subdomain(subdID));
                for (int i1 = 0; i1 < subdElementIds[subdID].GetLength(0); i1++)
                {
                    model.SubdomainsDictionary[subdID].Elements.Add(model.ElementsDictionary[subdElementIds[subdID][i1]]);//.ElementsDictionary.Add(subdElementIds[subdID][i1], model.ElementsDictionary[subdElementIds[subdID][i1]]);
                }
            }

        }

        public static void PrintDictionary (Dictionary<int, Dictionary<IDofType, int>> globalNodalDOFsDictionary,int TotalDOFs, int subdomainID)
        {
            double[] globalDOFs = new double[TotalDOFs];
            int counter = 0;

            foreach (int nodeID in globalNodalDOFsDictionary.Keys)
            {
                Dictionary<IDofType, int> dofTypes = globalNodalDOFsDictionary[nodeID];
                //Dictionary<DOFType, int> globalDOFTypes = new Dictionary<DOFType, int>(dofTypes.Count);
                foreach (IDofType dofType in dofTypes.Keys)
                {
                    if (dofTypes[dofType]!=-1)
                    {
                        globalDOFs[counter] = dofTypes[dofType];
                        counter += 1;
                    }
                }

            }

            string print_path_gen = @"C:\Users\turbo-x\Desktop\notes_elegxoi\MSOLVE_output_2\Subdomain{0}globalDOFs.txt";
            string file_no = subdomainID.ToString();
            string print_path = string.Format(print_path_gen, file_no);
            new Array1DWriter().WriteToFile(globalDOFs, print_path);
        }

        public static (Dictionary<int, Dictionary<int, IList<int>>>, Dictionary<int, List<int>>) FindEmbeddedElementsSubdomains(Model model)
        {
            //1
            Dictionary<int, Dictionary<int, IList<int>>> EmbeddedElementsHostSubdomainsAndSpecifcHexaElementsInThem = new Dictionary<int, Dictionary<int, IList<int>>>(); //embedded element- host subdomains -specific elements in subdomains
            //2
            Dictionary<int, List<int>> hexaConnectsShells = new Dictionary<int, List<int>>();
            //3
            List<int> totalEmbeddedElements = new List<int>();

            foreach (Element element in model.ElementsDictionary.Values)
            {
                if (element.ElementType is IEmbeddedElement)
                {
                    Dictionary<int, IList<int>> HostSubdomains = new Dictionary<int, IList<int>>();
                    foreach (var embeddedNode in ((IEmbeddedElement)element).EmbeddedNodes)
                    {
                        //1
                        Element hostELement = embeddedNode.EmbeddedInElement;
                        if (HostSubdomains.ContainsKey(hostELement.Subdomain.ID))
                        {
                            if (!HostSubdomains[hostELement.Subdomain.ID].Contains(hostELement.ID))
                            {
                                HostSubdomains[hostELement.Subdomain.ID].Add(hostELement.ID);
                            }
                        }
                        else
                        {
                            List<int> specificElementsIDs = new List<int>();
                            specificElementsIDs.Add(hostELement.ID);
                            HostSubdomains.Add(hostELement.Subdomain.ID, specificElementsIDs);                            
                        }
                        //2
                        if (hexaConnectsShells.ContainsKey(hostELement.ID))
                        {
                            if (!hexaConnectsShells[hostELement.ID].Contains(element.ID))
                            {
                                hexaConnectsShells[hostELement.ID].Add(element.ID);
                            }
                        }
                        else
                        {
                            List<int> connectionElementsData1 = new List<int>();
                            connectionElementsData1.Add(element.ID);
                            hexaConnectsShells.Add(hostELement.ID, connectionElementsData1);
                        }
                    }
                    if (HostSubdomains.Count > 1) // gia =1 den exoume dilhma gia to se poia subdomain tha entaxthei
                    { EmbeddedElementsHostSubdomainsAndSpecifcHexaElementsInThem.Add(element.ID, HostSubdomains); }
                }
            }
            return (EmbeddedElementsHostSubdomainsAndSpecifcHexaElementsInThem,hexaConnectsShells);
        }

        public static Dictionary<int, Dictionary<int, IList<int>>> FindAmbiguousEmbeddedElementsSubdomains(Model model)
        {
            Dictionary<int, Dictionary<int, IList<int>>> EmbeddedElementsHostSubdomainsAndElements = new Dictionary<int, Dictionary<int, IList<int>>>();

            foreach (Element element in model.ElementsDictionary.Values)
            {
                if (element.ElementType is IEmbeddedElement)
                {
                    Dictionary<int, IList<int>> HostSubdomains = new Dictionary<int, IList<int>>();
                    foreach (var embeddedNode in ((IEmbeddedElement)element).EmbeddedNodes)
                    {
                        Element hostELement = embeddedNode.EmbeddedInElement;
                        if (HostSubdomains.ContainsKey(hostELement.Subdomain.ID))
                        {
                            if (!HostSubdomains[hostELement.Subdomain.ID].Contains(hostELement.ID))
                            {
                                HostSubdomains[hostELement.Subdomain.ID].Add(hostELement.ID);
                            }
                        }
                        else
                        {
                            List<int> specificElementsIDs = new List<int>();
                            specificElementsIDs.Add(hostELement.ID);
                            HostSubdomains.Add(hostELement.Subdomain.ID, specificElementsIDs);
                        }
                    }

                    if (HostSubdomains.Count > 1)
                    { EmbeddedElementsHostSubdomainsAndElements.Add(element.ID, HostSubdomains); }
                }
            }
            return EmbeddedElementsHostSubdomainsAndElements;
        }

        public static Dictionary<int, List<int>> DetermineOnlyNeededCombinations( 
            Dictionary<int, Dictionary<int, IList<int>>> EmbeddedElementsHostSubdomainsAndSpecifcHexaElementsInThem,
            Dictionary<int, List<int>> hexaConnectsShells)
        {
            //List<List<>> CombinationElements
            Dictionary<int, List<int>> connectedShellElementsLists = new Dictionary<int, List<int>>();

            //List<int> connectedShells1 = new List<int>();

            foreach (int hexaID in hexaConnectsShells.Keys)
            {
                //foreach (int connectedShell in hexaConnectsShells[hexaID])
                //{
                List<int> foundInLists = new List<int>();

                foreach (int ListID in connectedShellElementsLists.Keys)
                {
                    bool isShellFoundInList = false;
                    
                    foreach (int shellELement in hexaConnectsShells[hexaID])
                    {
                        if (connectedShellElementsLists[ListID].Contains(shellELement))
                        {
                            isShellFoundInList = true;

                            if (!foundInLists.Contains(ListID))
                            { foundInLists.Add(ListID); }
                            //prosthese kai ta upoloipa shell sth lista
                            //kai oles tis upoloipes listes pou ta periexoun 
                        }
                    }
                }

                int foundInListsNum = foundInLists.Count();
                if (foundInListsNum==0)
                {
                    List<int> newConnectedShellsList = new List<int>();
                    foreach (int shellELement in hexaConnectsShells[hexaID])
                    {
                        newConnectedShellsList.Add(shellELement);
                    }

                    connectedShellElementsLists.Add(connectedShellElementsLists.Count()+1, newConnectedShellsList);
                }
                if (foundInListsNum==1)
                {
                    var updatedList = connectedShellElementsLists[foundInLists.ElementAt(0)];
                    foreach (int shellELement in hexaConnectsShells[hexaID])
                    {
                        if(!updatedList.Contains(shellELement))
                        { updatedList.Add(shellELement); }

                    }
                }
                if (foundInListsNum>1)
                {
                    //lists[1] = lists[1].Union(lists[2]).ToList();
                    //lists.Remove(2);

                    for (int i1=1; i1<foundInListsNum; i1++)
                    {
                        connectedShellElementsLists[foundInLists.ElementAt(0)] = 
                            connectedShellElementsLists[foundInLists.ElementAt(0)].Union(connectedShellElementsLists[foundInLists.ElementAt(i1)]).ToList();
                        //todo: concat can be used as well if it is known that there are not duplicates
                        connectedShellElementsLists.Remove(foundInLists.ElementAt(i1));
                    }
                }



                //connectedShellElementsLists[0].Union(connectedShellElementsLists[1]);

                //}
            }

            return connectedShellElementsLists;

        }

        public static void CalculateCombinationSolution(List<int> connectedShellElementsLists, Dictionary<int, Dictionary<int, IList<int>>> EmbeddedElementsHostSubdomainsAndElements)
        {
            int solutionVectorSize = connectedShellElementsLists.Count();
            int possibleSolutions = 1;            
            foreach(int shellId in connectedShellElementsLists)
            {
                possibleSolutions *= EmbeddedElementsHostSubdomainsAndElements[shellId].Count();
            }
            List<int>[] possibleSolutionHexas = new List<int>[possibleSolutions];
            //element ids twn hexas pou tha einai boundary

            int previousDivider = 1;
            int choices = 2;

            foreach (int shellId in connectedShellElementsLists)
            {
                for (int i1 = 0; i1 < solutionVectorSize; i1++)
                {
                    int SubregionsSize = solutionVectorSize / previousDivider;
                    int subregionsSeparationSize = SubregionsSize / EmbeddedElementsHostSubdomainsAndElements[shellId].Count();
                }
            }

            //var solutionhexaElementsIds = 
            int[] hexaIds = new int[2] { 0, 1 };
            var unique = hexaIds.Distinct().Count();

            int maxSubdElements = 0;
            foreach (var subdmNhexas in EmbeddedElementsHostSubdomainsAndElements.Values)
            {
                foreach(var hexaList in subdmNhexas.Values)
                {

                }
            }

        }

        public static int[] CalculateCombinationSolution2(List<int> connectedShellElementsList, Dictionary<int, Dictionary<int, IList<int>>> EmbeddedElementsHostSubdomainsAndElements)
        {
            int solutionVectorSize = connectedShellElementsList.Count();
            int possibleSolutions = 1;
            foreach (int shellId in connectedShellElementsList)
            {
                possibleSolutions *= EmbeddedElementsHostSubdomainsAndElements[shellId].Count();
            }
            List<int>[] possibleSolutionPenaltyHexas = new List<int>[possibleSolutions];
            //element ids twn hexas pou tha einai boundary




            int previousDivider = 1;

            foreach (int shellId in connectedShellElementsList)
            {
                //create chosenSubdomainPenaltyHexas : foreach subdomainId where the shellId will be added it contains the (host) hexaIds
                // that belong to the other not chosen subdomains
                Dictionary<int, List<int>> chosenSubdomainPenaltyHexas = new Dictionary<int, List<int>>();
                foreach (int ChosenSubdomainId in EmbeddedElementsHostSubdomainsAndElements[shellId].Keys)
                {
                    List<int> chosenSubdomainPenaltyHexaIDs = new List<int>();
                    foreach (int otherSubdomainsId in EmbeddedElementsHostSubdomainsAndElements[shellId].Keys)
                    {
                        if (!(otherSubdomainsId == ChosenSubdomainId))
                        {
                            foreach (int hexaID in EmbeddedElementsHostSubdomainsAndElements[shellId][otherSubdomainsId])
                            {
                                //edw thelei kai if hexaID den einai boundary logw allou shellId pou den einai ambiguous
                                if (!chosenSubdomainPenaltyHexaIDs.Contains(hexaID)) chosenSubdomainPenaltyHexaIDs.Add(hexaID);
                            }
                        }
                    }
                    chosenSubdomainPenaltyHexas.Add(ChosenSubdomainId, chosenSubdomainPenaltyHexaIDs);
                }

                // modify solution penalty hexas in possibleSolutionPenaltyHexas
                int regionsSize = possibleSolutions / previousDivider;
                int addedDivider = EmbeddedElementsHostSubdomainsAndElements[shellId].Count();
                int subregionsSize = regionsSize / addedDivider;


                for (int i2 = 0; i2 < previousDivider; i2++)
                {
                    for (int i3 = 0; i3 < addedDivider; i3++)
                    {
                        int subregionPosition = i2 * regionsSize + i3 * subregionsSize;
                        for (int i4 = 0; i4 < subregionsSize; i4++)
                        {
                            int solutionId = subregionPosition + i4;
                            if (possibleSolutionPenaltyHexas[solutionId] == null)
                            {
                                possibleSolutionPenaltyHexas[solutionId] = new List<int>();
                            }
                            possibleSolutionPenaltyHexas[solutionId] = possibleSolutionPenaltyHexas[solutionId].Union(chosenSubdomainPenaltyHexas.ElementAt(i3).Value).ToList();
                        }
                    }
                }
                previousDivider = previousDivider * addedDivider;


            }

            //count penalty hexas and estimate costs
            var possibleSolutionCosts = new int[possibleSolutions];
            for (int i1 = 0; i1 < possibleSolutions; i1++)
            {
                possibleSolutionCosts[i1] = possibleSolutionPenaltyHexas[i1].Count();
            }
            int chosenSolutionValue = possibleSolutionCosts.Min();
            int solutionPosition = 0;
            for (int i1 = 0; i1 < possibleSolutions; i1++)
            {
                if (possibleSolutionCosts[i1] == chosenSolutionValue)
                {
                    solutionPosition = i1;
                    break;
                }
            }

            //retrieve solution vector (subdomain IDs)
            var solutionVectorSubdomainIDs = new int[solutionVectorSize];
            previousDivider = 1;
            int remainder = solutionPosition;
            int counter = 0;
            int divider = 1;
            foreach (int shellId in connectedShellElementsList)
            {
                divider = divider * EmbeddedElementsHostSubdomainsAndElements[shellId].Count();
                int currentRegionSize = possibleSolutions / divider;
                int SolutionValuePosition = remainder / currentRegionSize;
                remainder = remainder % currentRegionSize;
                solutionVectorSubdomainIDs[counter] = EmbeddedElementsHostSubdomainsAndElements[shellId].ElementAt(SolutionValuePosition).Key;
                counter += 1;
            }

            return solutionVectorSubdomainIDs;
        }
    }

    /// <summary>
    /// Model separation methods for models with embedded elements
    /// Authors: Gerasimos Sotiropoulos
    /// </summary>
    public static class DdmCalculationsPartb // exei allaxei mono to onoma den exei ginei update se model xrhsh klp.
    {
        public static (Dictionary<int, Dictionary<int, IList<int>>>, Dictionary<int, List<int>>, Dictionary<int, List<int>>) FindEmbeddedElementsSubdomains(Model model, int totalSubdomains)
        {
            Dictionary<int, List<int>> AssignedSubdomains = new Dictionary<int, List<int>>(totalSubdomains);//TODO mporoume na tou dwsoume arxikh diastash ean thn exoume
            // to exw int (tou Dict dld) sumvolizei to subdomain ID
            // ta mesa int (dld afta pou periexei to List) einai ta IDs twn element pou tha mpoun se afth th subdomain

            //1
            Dictionary<int, Dictionary<int, IList<int>>> AmbiguousEmbeddedElementsHostSubdomainsAndSpecifcHexaElementsInThem = new Dictionary<int, Dictionary<int, IList<int>>>(); //embedded element- host subdomains -specific elements in subdomains
            // einai ola ta ambiguous

            //2
            Dictionary<int, List<int>> hexaConnectsShells = new Dictionary<int, List<int>>();
            //3
            List<int> totalEmbeddedElements = new List<int>();

            foreach (Element element in model.ElementsDictionary.Values) // ean xeroume apo thn arxh to id ton embedded mporoume na ton dinoume
            {
                if (element.ElementType is IEmbeddedElement)
                {
                    var e1 = element.ElementType as IEmbeddedElement;
                    Dictionary<int, IList<int>> HostSubdomains = new Dictionary<int, IList<int>>();
                    foreach (var embeddedNode in (e1).EmbeddedNodes)
                    {
                        //1
                        Element hostELement = embeddedNode.EmbeddedInElement;
                        if (HostSubdomains.ContainsKey(hostELement.Subdomain.ID))
                        {
                            if (!HostSubdomains[hostELement.Subdomain.ID].Contains(hostELement.ID))
                            {
                                HostSubdomains[hostELement.Subdomain.ID].Add(hostELement.ID);
                            }
                        }
                        else
                        {
                            List<int> specificElementsIDs = new List<int>();
                            specificElementsIDs.Add(hostELement.ID);
                            HostSubdomains.Add(hostELement.Subdomain.ID, specificElementsIDs);
                        }
                        //2
                        if (hexaConnectsShells.ContainsKey(hostELement.ID))
                        {
                            if (!hexaConnectsShells[hostELement.ID].Contains(element.ID))
                            {
                                hexaConnectsShells[hostELement.ID].Add(element.ID);
                            }
                        }
                        else
                        {
                            List<int> connectionElementsData1 = new List<int>();
                            connectionElementsData1.Add(element.ID);
                            hexaConnectsShells.Add(hostELement.ID, connectionElementsData1);
                        }
                    }
                    if (HostSubdomains.Count > 1) // gia =1 den exoume dilhma gia to se poia subdomain tha entaxthei
                    { AmbiguousEmbeddedElementsHostSubdomainsAndSpecifcHexaElementsInThem.Add(element.ID, HostSubdomains); }
                    if (HostSubdomains.Count == 1)
                    {
                        if (AssignedSubdomains.ContainsKey(HostSubdomains.ElementAt(0).Key))
                        {
                            AssignedSubdomains[HostSubdomains.ElementAt(0).Key].Add(element.ID);
                        }
                        else
                        {
                            List<int> subdElementsIds = new List<int>();
                            subdElementsIds.Add(element.ID);
                            AssignedSubdomains.Add(HostSubdomains.ElementAt(0).Key, subdElementsIds);
                        }
                    }

                }
            }
            return (AmbiguousEmbeddedElementsHostSubdomainsAndSpecifcHexaElementsInThem, hexaConnectsShells, AssignedSubdomains);
        }

        public static (Dictionary<int, Dictionary<int, IList<int>>>, Dictionary<int, List<int>>, Dictionary<int, List<int>>) FindEmbeddedElementsSubdomainsCorrected(Model model, int totalSubdomains)
        {
            Dictionary<int, List<int>> AssignedSubdomains = new Dictionary<int, List<int>>(totalSubdomains);//TODO mporoume na tou dwsoume arxikh diastash ean thn exoume
            // to exw int (tou Dict dld) sumvolizei to subdomain ID
            // ta mesa int (dld afta pou periexei to List) einai ta IDs twn element pou tha mpoun se afth th subdomain

            //1
            Dictionary<int, Dictionary<int, IList<int>>> AmbiguousEmbeddedElementsHostSubdomainsAndSpecifcHexaElementsInThem = new Dictionary<int, Dictionary<int, IList<int>>>(); //embedded element- host subdomains -specific elements in subdomains
            // einai ola ta ambiguous

            //2
            Dictionary<int, List<int>> hexaConnectsShells = new Dictionary<int, List<int>>();
            //3
            List<int> totalEmbeddedElements = new List<int>();

            foreach (Element element in model.ElementsDictionary.Values) // ean xeroume apo thn arxh to id ton embedded mporoume na ton dinoume
            {
                if (element.ElementType is IEmbeddedElement)
                {
                    Dictionary<int, List<int>> hexaConnectsShellsLocal = new Dictionary<int, List<int>>();

                    var e1 = element.ElementType as IEmbeddedElement;
                    Dictionary<int, IList<int>> HostSubdomains = new Dictionary<int, IList<int>>();
                    foreach (var embeddedNode in (e1).EmbeddedNodes)
                    {
                        //1
                        Element hostELement = embeddedNode.EmbeddedInElement;
                        if (HostSubdomains.ContainsKey(hostELement.Subdomain.ID))
                        {
                            if (!HostSubdomains[hostELement.Subdomain.ID].Contains(hostELement.ID))
                            {
                                HostSubdomains[hostELement.Subdomain.ID].Add(hostELement.ID);
                            }
                        }
                        else
                        {
                            List<int> specificElementsIDs = new List<int>();
                            specificElementsIDs.Add(hostELement.ID);
                            HostSubdomains.Add(hostELement.Subdomain.ID, specificElementsIDs);
                        }
                        //2
                        if (hexaConnectsShellsLocal.ContainsKey(hostELement.ID))
                        {
                            if (!hexaConnectsShellsLocal[hostELement.ID].Contains(element.ID))
                            {
                                hexaConnectsShellsLocal[hostELement.ID].Add(element.ID);
                            }
                        }
                        else
                        {
                            List<int> connectionElementsData1 = new List<int>();
                            connectionElementsData1.Add(element.ID);
                            hexaConnectsShellsLocal.Add(hostELement.ID, connectionElementsData1);
                        }
                    }
                    if (HostSubdomains.Count > 1) // gia =1 den exoume dilhma gia to se poia subdomain tha entaxthei
                    {
                        AmbiguousEmbeddedElementsHostSubdomainsAndSpecifcHexaElementsInThem.Add(element.ID, HostSubdomains);
                        foreach (int hexaID in hexaConnectsShellsLocal.Keys)
                        {
                            if (hexaConnectsShells.Keys.Contains(hexaID))
                            {
                                foreach (int connectedShellId in hexaConnectsShellsLocal[hexaID])
                                {
                                    if (!hexaConnectsShells[hexaID].Contains(connectedShellId))
                                    {
                                        hexaConnectsShells[hexaID].Add(connectedShellId);
                                    }
                                }
                            }
                            else
                            {
                                hexaConnectsShells.Add(hexaID, new List<int>());
                                foreach (int connectedShellId in hexaConnectsShellsLocal[hexaID])
                                {
                                    hexaConnectsShells[hexaID].Add(connectedShellId);
                                }
                            }
                        }
                    }
                    if (HostSubdomains.Count == 1)
                    {
                        if (AssignedSubdomains.ContainsKey(HostSubdomains.ElementAt(0).Key))
                        {
                            AssignedSubdomains[HostSubdomains.ElementAt(0).Key].Add(element.ID);
                        }
                        else
                        {
                            List<int> subdElementsIds = new List<int>();
                            subdElementsIds.Add(element.ID);
                            AssignedSubdomains.Add(HostSubdomains.ElementAt(0).Key, subdElementsIds);
                        }
                    }

                }
            }
            return (AmbiguousEmbeddedElementsHostSubdomainsAndSpecifcHexaElementsInThem, hexaConnectsShells, AssignedSubdomains);
        }

        public static int[][] DetermineAmbiguousSimple(Dictionary<int, Dictionary<int, IList<int>>> AmbiguousEmbeddedElementsHostSubdomainsAndSpecifcHexaElementsInThem, Dictionary<int, List<int>> hexaConnectsShells, Dictionary<int, List<int>> AssignedSubdomains)
        {
            Dictionary<int, Dictionary<int, IList<int>>> AmbigElements = AmbiguousEmbeddedElementsHostSubdomainsAndSpecifcHexaElementsInThem;

            foreach (int ambElementID in AmbigElements.Keys)
            {
                int[] numSpecificElements = new int[AmbigElements[ambElementID].Keys.Count];

                int subdCounter = 0;
                foreach (int subdID in AmbigElements[ambElementID].Keys)
                {
                    //numSpecificElements[subdCounter] = AmbigElements[ambElementID][subdID].Count; // Count dld ta element ths inner listas tou AmbigElements...
                    numSpecificElements[subdCounter] = AmbigElements[ambElementID].ElementAt(subdCounter).Value.Count; // Value einai h inner lista tou AmbigElements...
                    subdCounter++;
                }

                int thesiChosenSubd = 0;
                int nElements = 0;
                for (int i1 = 0; i1 < numSpecificElements.GetLength(0); i1++)
                {
                    if (numSpecificElements[i1] > nElements)
                    {
                        thesiChosenSubd = i1;
                        nElements = numSpecificElements[i1];
                    }
                }

                int chosenSubdomainID = AmbigElements[ambElementID].ElementAt(thesiChosenSubd).Key;

                if (AssignedSubdomains.Keys.Contains(chosenSubdomainID))
                {
                    AssignedSubdomains[chosenSubdomainID].Add(ambElementID);
                }
                else
                {
                    List<int> subdElementsIds = new List<int>();
                    subdElementsIds.Add(ambElementID);
                    AssignedSubdomains.Add(chosenSubdomainID, subdElementsIds);
                }

            }

            return ConvertIntListToArray(AssignedSubdomains);
        }

        public static int[][] ConvertIntListToArray(Dictionary<int, List<int>> AssignedSubdomains)
        {
            int maxSubdId = AssignedSubdomains.Keys.Max();
            int[][] subdIdsAndElements = new int[maxSubdId + 1][]; //todo:revisit this

            foreach (int subdID in AssignedSubdomains.Keys)
            {
                subdIdsAndElements[subdID] = AssignedSubdomains[subdID].ToArray();
            }

            return subdIdsAndElements;
        }

        //public static void MakeModelDictionariesZeroBasedForDecomposer(Model model)
        //{
        //    model.SubdomainsDictionary[1].ID = 0;
        //    Subdomain subdomain_ini = model.SubdomainsDictionary[1];
        //    model.SubdomainsDictionary.Remove(1);
        //    model.SubdomainsDictionary.Add(0, subdomain_ini);

        //    Dictionary<int, Element> ElementsDictionary_2 = new Dictionary<int, Element>(model.ElementsDictionary.Count);
        //    for (int i1 = 0; i1 < model.ElementsDictionary.Count; i1++)
        //    {
        //        ElementsDictionary_2.Add(model.ElementsDictionary[i1 + 1].ID - 1, model.ElementsDictionary[i1 + 1]);
        //        ElementsDictionary_2[model.ElementsDictionary[i1 + 1].ID - 1].ID += -1;
        //    }
        //    int nElement = model.ElementsDictionary.Count;
        //    for (int i1 = 0; i1 < nElement; i1++)
        //    {
        //        model.ElementsDictionary.Remove(i1 + 1);
        //        model.SubdomainsDictionary[0].ElementsDictionary.Remove(i1 + 1);
        //    }
        //    for (int i1 = 0; i1 < nElement; i1++)
        //    {
        //        model.ElementsDictionary.Add(ElementsDictionary_2[i1].ID, ElementsDictionary_2[i1]);
        //        model.SubdomainsDictionary[0].ElementsDictionary.Add(ElementsDictionary_2[i1].ID, ElementsDictionary_2[i1]);
        //    }

        //    Dictionary<int, Node> NodesDictionary_2 = new Dictionary<int, Node>(model.NodesDictionary.Count);
        //    for (int i1 = 0; i1 < model.NodesDictionary.Count; i1++)
        //    {
        //        NodesDictionary_2.Add(model.NodesDictionary[i1 + 1].ID - 1, model.NodesDictionary[i1 + 1]);
        //        NodesDictionary_2[model.NodesDictionary[i1 + 1].ID - 1].ID += -1;
        //    }

        //    int nNode = model.NodesDictionary.Count;
        //    for (int i1 = 0; i1 < nNode; i1++)
        //    {
        //        model.NodesDictionary.Remove(i1 + 1);
        //        //model.SubdomainsDictionary[0].NodesDictionary.Remove(i1 + 1); // den peirazoume to subdomain nodesDictionary, ftiahnetai mono tou (pithanws apo to connect data Structures)
        //    }
        //    for (int i1 = 0; i1 < nNode; i1++)
        //    {
        //        model.NodesDictionary.Add(NodesDictionary_2[i1].ID, NodesDictionary_2[i1]);
        //        //model.SubdomainsDictionary[0].NodesDictionary.Add(NodesDictionary_2[i1].ID, NodesDictionary_2[i1]); // den peirazoume to subdomain nodesDictionary, ftiahnetai mono tou (pithanws apo to connect data Structures)
        //    }
        //}

        public static Dictionary<int, List<int>> FindEmbeddedElementsSubdomainsCorrectedSimple(Model model, int totalSubdomains)
        {
            Dictionary<int, List<int>> AssignedSubdomains = new Dictionary<int, List<int>>(totalSubdomains);//TODO mporoume na tou dwsoume arxikh diastash ean thn exoume
            // to exw int (tou Dict dld) sumvolizei to subdomain ID
            // ta mesa int (dld afta pou periexei to List) einai ta IDs twn element pou tha mpoun se afth th subdomain

            //1
            //Dictionary<int, Dictionary<int, IList<int>>> AmbiguousEmbeddedElementsHostSubdomainsAndSpecifcHexaElementsInThem = new Dictionary<int, Dictionary<int, IList<int>>>(); //embedded element- host subdomains -specific elements in subdomains
            // einai ola ta ambiguous

            //2
            //Dictionary<int, List<int>> hexaConnectsShells = new Dictionary<int, List<int>>();
            //3
            List<int> totalEmbeddedElements = new List<int>();

            foreach (Element element in model.ElementsDictionary.Values) // ean xeroume apo thn arxh to id ton embedded mporoume na ton dinoume
            {
                if (element.ElementType is IEmbeddedElement)
                {
                    Dictionary<int, List<int>> hexaConnectsShellsLocal = new Dictionary<int, List<int>>();

                    var e1 = element.ElementType as IEmbeddedElement;
                    Dictionary<int, IList<int>> HostSubdomains = new Dictionary<int, IList<int>>();
                    foreach (var embeddedNode in (e1).EmbeddedNodes)
                    {
                        //1
                        Element hostELement = embeddedNode.EmbeddedInElement;
                        if (HostSubdomains.ContainsKey(hostELement.Subdomain.ID))
                        {
                            if (!HostSubdomains[hostELement.Subdomain.ID].Contains(hostELement.ID))
                            {
                                HostSubdomains[hostELement.Subdomain.ID].Add(hostELement.ID);
                            }
                        }
                        else
                        {
                            List<int> specificElementsIDs = new List<int>();
                            specificElementsIDs.Add(hostELement.ID);
                            HostSubdomains.Add(hostELement.Subdomain.ID, specificElementsIDs);
                        }
                        //2
                        //if (hexaConnectsShellsLocal.ContainsKey(hostELement.ID))
                        //{
                        //    if (!hexaConnectsShellsLocal[hostELement.ID].Contains(element.ID))
                        //    {
                        //        hexaConnectsShellsLocal[hostELement.ID].Add(element.ID);
                        //    }
                        //}
                        //else
                        //{
                        //    List<int> connectionElementsData1 = new List<int>();
                        //    connectionElementsData1.Add(element.ID);
                        //    hexaConnectsShellsLocal.Add(hostELement.ID, connectionElementsData1);
                        //}
                    }
                    if (HostSubdomains.Count > 1) // gia =1 den exoume dilhma gia to se poia subdomain tha entaxthei
                    {
                        int chosenSubdomainId = 0;
                        int hexaListlength = 0;
                        foreach (int subdId in HostSubdomains.Keys)
                        {
                            if (HostSubdomains[subdId].Count > hexaListlength)
                            {
                                chosenSubdomainId = subdId;
                                hexaListlength = HostSubdomains[subdId].Count;
                            }
                        }
                        if (AssignedSubdomains.ContainsKey(chosenSubdomainId))
                        {
                            AssignedSubdomains[chosenSubdomainId].Add(element.ID);
                        }
                        else
                        {
                            List<int> subdElementsIds = new List<int>();
                            subdElementsIds.Add(element.ID);
                            AssignedSubdomains.Add(chosenSubdomainId, subdElementsIds);
                        }
                    }
                    if (HostSubdomains.Count == 1)
                    {
                        if (AssignedSubdomains.ContainsKey(HostSubdomains.ElementAt(0).Key))
                        {
                            AssignedSubdomains[HostSubdomains.ElementAt(0).Key].Add(element.ID);
                        }
                        else
                        {
                            List<int> subdElementsIds = new List<int>();
                            subdElementsIds.Add(element.ID);
                            AssignedSubdomains.Add(HostSubdomains.ElementAt(0).Key, subdElementsIds);
                        }
                    }

                }
            }
            return AssignedSubdomains;
        }
    }

    /// <summary>
    /// Model separation methods for models with embedded elements
    /// Authors: Gerasimos Sotiropoulos
    /// </summary>
    public static class DdmCalculationsGeneral
    {
        public static void BuildModelInterconnectionData(Model model) //xreiazetai na ginei oloklhro update se model xrhsh dioti exei prosarmostei mono h methodos h prwth
        {
            //private void BuildSubdomainOfEachElement()
            foreach (Subdomain subdomain in model.SubdomainsDictionary.Values)
            {
                foreach (Element element in subdomain.Elements)
                { element.Subdomain = subdomain; }
            }

            //private void BuildElementDictionaryOfEachNode()            
            foreach (Element element in model.ElementsDictionary.Values)
            {
                foreach (Node node in element.Nodes)
                { node.ElementsDictionary.Add(element.ID, element); }
            }

            foreach (Node node in model.NodesDictionary.Values)
            { node.BuildSubdomainDictionary(); }


            //TEMP COMMENT OUT
            //foreach (Subdomain subdomain in model.SubdomainsDictionary.Values)
            //{ subdomain.DefineNodesFromElements(); }  //BuildNodesDictionary(); }
            ////TODO afto tha ginei commented out afou den mporoume na to kanoume undo build meta
        }

        public static void UndoModelInterconnectionDataBuild(Model model)
        {
            //private void BuildSubdomainOfEachElement()
            foreach (Subdomain subdomain in model.SubdomainsDictionary.Values)
            {
                foreach (Element element in subdomain.Elements)
                { element.Subdomain = null; }  // subdomain; }
            }

            //private void BuildElementDictionaryOfEachNode()            
            //foreach (Element element in model.ElementsDictionary.Values)
            //{
            //    foreach (Node node in element.Nodes)
            //    { node.ElementsDictionary.Add(element.ID, element); }
            //}
            foreach (Node node in model.NodesDictionary.Values)
            {
                node.ElementsDictionary.Clear();
                node.SubdomainsDictionary.Clear(); //to ena mono egine v2 opws fainetai sto Node.BuildSubdomainDictionary().
            }
            //foreach (Node node in model.NodesDictionary.Values)
            //{ node.BuildSubdomainDictionary(); }


            //TEMP COMMENT OUT opote den xreiazetai na ginei update se v2
            //foreach (Subdomain subdomain in model.SubdomainsDictionary.Values)
            //{
            //    //entolh antistoixh tou exhs: //subdomain.BuildNodesDictionary();
            //    subdomain.NodesDictionary.Clear();
            //}
        }

        public static int[][] DetermineHexaElementsSubdomainsFromModel(Model model)
        {
            int[][] subdomainsAndHexas = new int[model.Subdomains.Count()][];

            for (int subdomainId = 0; subdomainId < model.Subdomains.Count(); subdomainId++)
            {
                var subdomain = model.Subdomains[subdomainId]; //ZERo based model.subdomainsDictionary access == model.Subdomains access
                subdomainsAndHexas[subdomainId] = new int[subdomain.Elements.Count()];
                int hexaPositionInArray = 0;
                foreach (Element element in subdomain.Elements)
                {
                    subdomainsAndHexas[subdomainId][hexaPositionInArray] = element.ID;
                    hexaPositionInArray++;
                }
            }

            return subdomainsAndHexas;
        }

        public static int[][] CombineSubdomainElementsIdsArraysIntoOne(int[][] subdomainsAndElements1, int[][] subdomainsAndElements2)
        {
            // input example: int[][] subdomainsAndHexas, int[][] subdomainsAndEmbedded
            int[][] subdomainsAndElements = new int[Math.Max(subdomainsAndElements1.GetLength(0), subdomainsAndElements2.GetLength(0))][];

            if (subdomainsAndElements1.GetLength(0) == subdomainsAndElements2.GetLength(0))
            {
                for (int i1 = 0; i1 < subdomainsAndElements1.GetLength(0); i1++)
                {
                    if (subdomainsAndElements1[i1] == null)
                    {
                        if (!(subdomainsAndElements2[i1] == null))
                        {
                            subdomainsAndElements[i1] = new int[subdomainsAndElements2[i1].Length];
                            for (int i2 = 0; i2 < subdomainsAndElements2[i1].Length; i2++) { subdomainsAndElements[i1][i2] = subdomainsAndElements2[i1][i2]; }
                        }
                    }
                    if (subdomainsAndElements2[i1] == null)
                    {
                        if (!(subdomainsAndElements1[i1] == null))
                        {
                            subdomainsAndElements[i1] = new int[subdomainsAndElements1[i1].Length];
                            for (int i2 = 0; i2 < subdomainsAndElements1[i1].Length; i2++) { subdomainsAndElements[i1][i2] = subdomainsAndElements1[i1][i2]; }
                        }
                    }
                    if (!(subdomainsAndElements1[i1] == null))
                    {
                        if (!(subdomainsAndElements2[i1] == null))
                        {
                            subdomainsAndElements[i1] = subdomainsAndElements1[i1].Union(subdomainsAndElements2[i1]).ToArray();
                        }
                    }

                }
            }

            if (subdomainsAndElements1.GetLength(0) < subdomainsAndElements2.GetLength(0))
            {
                for (int i1 = 0; i1 < subdomainsAndElements1.GetLength(0); i1++)
                {
                    if (subdomainsAndElements1[i1] == null)
                    {
                        if (!(subdomainsAndElements2[i1] == null))
                        {
                            subdomainsAndElements[i1] = new int[subdomainsAndElements2[i1].Length];
                            for (int i2 = 0; i2 < subdomainsAndElements2[i1].Length; i2++) { subdomainsAndElements[i1][i2] = subdomainsAndElements2[i1][i2]; }
                        }
                    }
                    if (subdomainsAndElements2[i1] == null)
                    {
                        if (!(subdomainsAndElements1[i1] == null))
                        {
                            subdomainsAndElements[i1] = new int[subdomainsAndElements1[i1].Length];
                            for (int i2 = 0; i2 < subdomainsAndElements1[i1].Length; i2++) { subdomainsAndElements[i1][i2] = subdomainsAndElements1[i1][i2]; }
                        }
                    }
                    if (!(subdomainsAndElements1[i1] == null))
                    {
                        if (!(subdomainsAndElements2[i1] == null))
                        {
                            subdomainsAndElements[i1] = subdomainsAndElements1[i1].Union(subdomainsAndElements2[i1]).ToArray();
                        }
                    }

                }

                for (int i1 = subdomainsAndElements1.GetLength(0); i1 < subdomainsAndElements2.GetLength(0); i1++)
                {
                    subdomainsAndElements[i1] = new int[subdomainsAndElements2[i1].Length];
                    for (int i2 = 0; i2 < subdomainsAndElements2[i1].Length; i2++) { subdomainsAndElements[i1][i2] = subdomainsAndElements2[i1][i2]; }
                }
            }

            if (subdomainsAndElements1.GetLength(0) < subdomainsAndElements2.GetLength(0))
            {
                for (int i1 = 0; i1 < subdomainsAndElements2.GetLength(0); i1++)
                {
                    if (subdomainsAndElements1[i1] == null)
                    {
                        if (!(subdomainsAndElements2[i1] == null))
                        {
                            subdomainsAndElements[i1] = new int[subdomainsAndElements2[i1].Length];
                            for (int i2 = 0; i2 < subdomainsAndElements2[i1].Length; i2++) { subdomainsAndElements[i1][i2] = subdomainsAndElements2[i1][i2]; }
                        }
                    }
                    if (subdomainsAndElements2[i1] == null)
                    {
                        if (!(subdomainsAndElements1[i1] == null))
                        {
                            subdomainsAndElements[i1] = new int[subdomainsAndElements1[i1].Length];
                            for (int i2 = 0; i2 < subdomainsAndElements1[i1].Length; i2++) { subdomainsAndElements[i1][i2] = subdomainsAndElements1[i1][i2]; }
                        }
                    }
                    if (!(subdomainsAndElements1[i1] == null))
                    {
                        if (!(subdomainsAndElements2[i1] == null))
                        {
                            subdomainsAndElements[i1] = subdomainsAndElements1[i1].Union(subdomainsAndElements2[i1]).ToArray();
                        }
                    }

                }

                for (int i1 = subdomainsAndElements2.GetLength(0); i1 < subdomainsAndElements1.GetLength(0); i1++)
                {
                    subdomainsAndElements[i1] = new int[subdomainsAndElements1[i1].Length];
                    for (int i2 = 0; i2 < subdomainsAndElements1[i1].Length; i2++) { subdomainsAndElements[i1][i2] = subdomainsAndElements1[i1][i2]; }
                }
            }


            return subdomainsAndElements;
        }

        public static int[][] CombineSubdomainElementsIdsArraysIntoOneCopy(int[][] subdomainsAndElements1, int[][] subdomainsAndElements2)
        {
            // input example: int[][] subdomainsAndHexas, int[][] subdomainsAndEmbedded
            int[][] subdomainsAndElements = new int[subdomainsAndElements1.GetLength(0)][];

            for (int i1 = 0; i1 < subdomainsAndElements1.GetLength(0); i1++)
            {
                if (subdomainsAndElements1[i1] == null)
                {
                    if (!(subdomainsAndElements2[i1] == null))
                    {
                        subdomainsAndElements[i1] = new int[subdomainsAndElements2[i1].Length];
                        for (int i2 = 0; i2 < subdomainsAndElements2[i1].Length; i2++) { subdomainsAndElements[i1][i2] = subdomainsAndElements2[i1][i2]; }
                    }
                }
                if (subdomainsAndElements2[i1] == null)
                {
                    if (!(subdomainsAndElements1[i1] == null))
                    {
                        subdomainsAndElements[i1] = new int[subdomainsAndElements1[i1].Length];
                        for (int i2 = 0; i2 < subdomainsAndElements1[i1].Length; i2++) { subdomainsAndElements[i1][i2] = subdomainsAndElements1[i1][i2]; }
                    }
                }
                if (!(subdomainsAndElements1[i1] == null))
                {
                    if (!(subdomainsAndElements2[i1] == null))
                    {
                        subdomainsAndElements[i1] = subdomainsAndElements1[i1].Union(subdomainsAndElements2[i1]).ToArray();
                    }
                }

            }
            return subdomainsAndElements;
        }

        public static int[][] DetermineCoheiveELementsSubdomains(Model model, int totalSubdomains)
        {
            (Dictionary<int, Dictionary<int, IList<int>>> AmbiguousEmbeddedElementsHostSubdomainsAndSpecifcHexaElementsInThem,
               Dictionary<int, List<int>> hexaConnectsShells,
               Dictionary<int, List<int>> AssignedSubdomains) =
               DdmCalculationsPartb.FindEmbeddedElementsSubdomainsCorrected(model, totalSubdomains);

            Dictionary<int, List<int>> connectedShellElementsLists =
                DdmCalculations.DetermineOnlyNeededCombinations(
            AmbiguousEmbeddedElementsHostSubdomainsAndSpecifcHexaElementsInThem,
             hexaConnectsShells);

            int numlistsForCombinations = connectedShellElementsLists.Keys.Count();
            int[][] combinationSolutions = new int[numlistsForCombinations][];
            for (int i1 = 0; i1 < numlistsForCombinations; i1++)
            {
                combinationSolutions[i1] =
                    DdmCalculations.CalculateCombinationSolution2(connectedShellElementsLists.ElementAt(i1).Value,
                    AmbiguousEmbeddedElementsHostSubdomainsAndSpecifcHexaElementsInThem);
            }

            //gather solutions 
            var subdElementsNums = new int[totalSubdomains];
            for (int i1 = 0; i1 < numlistsForCombinations; i1++)
            {
                for (int i2 = 0; i2 < combinationSolutions[i1].Length; i2++)
                {
                    subdElementsNums[combinationSolutions[i1][i2]] += 1;
                }
            }

            int[][] subdAmbCohElementIds = new int[totalSubdomains][];
            for (int i1 = 0; i1 < totalSubdomains; i1++)
            {
                subdAmbCohElementIds[i1] = new int[subdElementsNums[i1]];
            }

            int[] subdCohElemCounters = new int[totalSubdomains];
            for (int i1 = 0; i1 < numlistsForCombinations; i1++)
            {
                var cohIDsList = connectedShellElementsLists.ElementAt(i1).Value;
                var solutionSubdomainIDs = combinationSolutions[i1];

                for (int i2 = 0; i2 < connectedShellElementsLists.ElementAt(i1).Value.Count; i2++)
                {
                    int shellID = cohIDsList[i2];
                    int subdId = solutionSubdomainIDs[i2];
                    subdAmbCohElementIds[subdId][subdCohElemCounters[subdId]] = shellID;
                    subdCohElemCounters[subdId] += 1;
                }
            }

            int[][] subdCohElementIdsDirect = DdmCalculationsPartb.ConvertIntListToArray(AssignedSubdomains);
            int[][] subdCohElementIds = CombineSubdomainElementsIdsArraysIntoOne(subdAmbCohElementIds, subdCohElementIdsDirect);

            return subdCohElementIds;
        }

        public static int[][] DetermineCoheiveELementsSubdomainsSimple(Model model, int totalSubdomains)
        {
            Dictionary<int, List<int>> AssignedSubdomains =
               DdmCalculationsPartb.FindEmbeddedElementsSubdomainsCorrectedSimple(model, totalSubdomains);

            int[][] subdCohElementIdsDirect = DdmCalculationsPartb.ConvertIntListToArray(AssignedSubdomains);
            return subdCohElementIdsDirect;
        }

        public static int[][] DetermineShellELementsSubdomains(Model model, int totalSubdomains, int[][] subdCohElementIds,
            int[] lowerCohesiveBound, int[] upperCohesiveBound, int[] grShElementssnumber)
        {
            List<int>[] subdShellElementIds = new List<int>[totalSubdomains];
            for (int i1 = 0; i1 < totalSubdomains; i1++)
            {
                if (!(subdCohElementIds[i1] == null))
                { subdShellElementIds[i1] = new List<int>(subdCohElementIds[i1].Length); }
            }

            for (int i1 = 0; i1 < totalSubdomains; i1++)
            {
                if (!(subdCohElementIds[i1] == null))
                {
                    for (int i2 = 0; i2 < subdCohElementIds[i1].Length; i2++)
                    {
                        int cohID = subdCohElementIds[i1][i2];
                        for (int i3 = 0; i3 < lowerCohesiveBound.Length; i3++)
                        {
                            if ((lowerCohesiveBound[i3] <= cohID) & (upperCohesiveBound[i3] >= cohID))
                            {
                                //subdID=i1;
                                if (!subdShellElementIds[i1].Contains(cohID - grShElementssnumber[i3]))
                                { subdShellElementIds[i1].Add(cohID - grShElementssnumber[i3]); }
                                break;
                            }

                        }
                    }
                }
            }

            int[][] subdShellElementIdsArrays = new int[subdShellElementIds.Length][];
            for (int i1 = 0; i1 < subdShellElementIds.Length; i1++)
            {
                if (!(subdShellElementIds[i1] == null))
                { subdShellElementIdsArrays[i1] = subdShellElementIds[i1].ToArray(); }
            }

            return subdShellElementIdsArrays;
        }

        public static void PrintSubdomainDataForPostPro(int[][] subdHexaIds, int[][] subdCohElementIds, int[][] subdShellElementIds, string generalPath)
        {
            string hexaPath = generalPath + @"\subdomainHexas.txt";
            string cohPath = generalPath + @"\subdomainCohesiveElements.txt";
            string shellPath = generalPath + @"\subdomainShellElements.txt";

            #region hexa elements
            int hexaPrintLength = 0;
            for (int i1 = 0; i1 < subdHexaIds.Length; i1++)
            {
                if (subdHexaIds[i1] == null)
                {
                    hexaPrintLength += 1;
                }
                else
                {
                    hexaPrintLength += 1 + subdHexaIds[i1].Length;
                }
            }

            var hexaPrint = new int[hexaPrintLength];
            int hexaPrintCounter = 0;
            for (int i1 = 0; i1 < subdHexaIds.Length; i1++)
            {
                if (subdHexaIds[i1] == null)
                {
                    hexaPrint[hexaPrintCounter] = 0;
                    hexaPrintCounter += 1;
                }
                else
                {
                    hexaPrint[hexaPrintCounter] = subdHexaIds[i1].Length;
                    hexaPrintCounter += 1;
                    for (int i2 = 0; i2 < subdHexaIds[i1].Length; i2++)
                    {
                        hexaPrint[hexaPrintCounter] = subdHexaIds[i1][i2];
                        hexaPrintCounter += 1;
                    }
                }
            }
            WriteToFileVector(hexaPrint, hexaPath);
            #endregion

            #region cohesive elements
            int cohePrintLength = 0;
            for (int i1 = 0; i1 < subdCohElementIds.Length; i1++)
            {
                if (subdCohElementIds[i1] == null)
                {
                    cohePrintLength += 1;
                }
                else
                {
                    cohePrintLength += 1 + subdCohElementIds[i1].Length;
                }
            }

            var cohePrint = new int[cohePrintLength];
            int cohePrintCounter = 0;
            for (int i1 = 0; i1 < subdCohElementIds.Length; i1++)
            {
                if (subdCohElementIds[i1] == null)
                {
                    cohePrint[cohePrintCounter] = 0;
                    cohePrintCounter += 1;
                }
                else
                {
                    cohePrint[cohePrintCounter] = subdCohElementIds[i1].Length;
                    cohePrintCounter += 1;
                    for (int i2 = 0; i2 < subdCohElementIds[i1].Length; i2++)
                    {
                        cohePrint[cohePrintCounter] = subdCohElementIds[i1][i2];
                        cohePrintCounter += 1;
                    }
                }
            }
            WriteToFileVector(cohePrint, cohPath);
            #endregion

            #region shell elements
            int shellPrintLength = 0;
            for (int i1 = 0; i1 < subdShellElementIds.Length; i1++)
            {
                if (subdShellElementIds[i1] == null)
                {
                    shellPrintLength += 1;
                }
                else
                {
                    shellPrintLength += 1 + subdShellElementIds[i1].Length;
                }
            }

            var shellPrint = new int[shellPrintLength];
            int shellPrintCounter = 0;
            for (int i1 = 0; i1 < subdShellElementIds.Length; i1++)
            {
                if (subdShellElementIds[i1] == null)
                {
                    shellPrint[shellPrintCounter] = 0;
                    shellPrintCounter += 1;
                }
                else
                {
                    shellPrint[shellPrintCounter] = subdShellElementIds[i1].Length;
                    shellPrintCounter += 1;
                    for (int i2 = 0; i2 < subdShellElementIds[i1].Length; i2++)
                    {
                        shellPrint[shellPrintCounter] = subdShellElementIds[i1][i2];
                        shellPrintCounter += 1;
                    }
                }
            }
            WriteToFileVector(shellPrint, shellPath);
            #endregion


        }

        public static (int[], int[], int[]) GetSubdomainDataForPostPro(int[][] subdHexaIds, int[][] subdCohElementIds, int[][] subdShellElementIds, string generalPath)
        {
            string hexaPath = generalPath + @"\subdomainHexas.txt";
            string cohPath = generalPath + @"\subdomainCohesiveElements.txt";
            string shellPath = generalPath + @"\subdomainShellElements.txt";

            #region hexa elements
            int hexaPrintLength = 0;
            for (int i1 = 0; i1 < subdHexaIds.Length; i1++)
            {
                if (subdHexaIds[i1] == null)
                {
                    hexaPrintLength += 1;
                }
                else
                {
                    hexaPrintLength += 1 + subdHexaIds[i1].Length;
                }
            }

            var hexaPrint = new int[hexaPrintLength];
            int hexaPrintCounter = 0;
            for (int i1 = 0; i1 < subdHexaIds.Length; i1++)
            {
                if (subdHexaIds[i1] == null)
                {
                    hexaPrint[hexaPrintCounter] = 0;
                    hexaPrintCounter += 1;
                }
                else
                {
                    hexaPrint[hexaPrintCounter] = subdHexaIds[i1].Length;
                    hexaPrintCounter += 1;
                    for (int i2 = 0; i2 < subdHexaIds[i1].Length; i2++)
                    {
                        hexaPrint[hexaPrintCounter] = subdHexaIds[i1][i2];
                        hexaPrintCounter += 1;
                    }
                }
            }
            //WriteToFileVector(hexaPrint, hexaPath);
            #endregion

            #region cohesive elements
            int cohePrintLength = 0;
            for (int i1 = 0; i1 < subdCohElementIds.Length; i1++)
            {
                if (subdCohElementIds[i1] == null)
                {
                    cohePrintLength += 1;
                }
                else
                {
                    cohePrintLength += 1 + subdCohElementIds[i1].Length;
                }
            }

            var cohePrint = new int[cohePrintLength];
            int cohePrintCounter = 0;
            for (int i1 = 0; i1 < subdCohElementIds.Length; i1++)
            {
                if (subdCohElementIds[i1] == null)
                {
                    cohePrint[cohePrintCounter] = 0;
                    cohePrintCounter += 1;
                }
                else
                {
                    cohePrint[cohePrintCounter] = subdCohElementIds[i1].Length;
                    cohePrintCounter += 1;
                    for (int i2 = 0; i2 < subdCohElementIds[i1].Length; i2++)
                    {
                        cohePrint[cohePrintCounter] = subdCohElementIds[i1][i2];
                        cohePrintCounter += 1;
                    }
                }
            }
            //WriteToFileVector(cohePrint, cohPath);
            #endregion

            #region shell elements
            int shellPrintLength = 0;
            for (int i1 = 0; i1 < subdShellElementIds.Length; i1++)
            {
                if (subdShellElementIds[i1] == null)
                {
                    shellPrintLength += 1;
                }
                else
                {
                    shellPrintLength += 1 + subdShellElementIds[i1].Length;
                }
            }

            var shellPrint = new int[shellPrintLength];
            int shellPrintCounter = 0;
            for (int i1 = 0; i1 < subdShellElementIds.Length; i1++)
            {
                if (subdShellElementIds[i1] == null)
                {
                    shellPrint[shellPrintCounter] = 0;
                    shellPrintCounter += 1;
                }
                else
                {
                    shellPrint[shellPrintCounter] = subdShellElementIds[i1].Length;
                    shellPrintCounter += 1;
                    for (int i2 = 0; i2 < subdShellElementIds[i1].Length; i2++)
                    {
                        shellPrint[shellPrintCounter] = subdShellElementIds[i1][i2];
                        shellPrintCounter += 1;
                    }
                }
            }
            //WriteToFileVector(shellPrint, shellPath);
            #endregion

            return (hexaPrint, cohePrint, shellPrint);
        }

        public static void WriteToFileVector(int[] array, string path2)
        {
            var writer2 = new StreamWriter(path2);
            for (int i = 0; i < array.GetLength(0); ++i)
            {
                writer2.Write(array[i]);
                writer2.Write(' ');
                writer2.WriteLine(); // allagh seiras (dld grafei oti exei mesa h parenths=esh edw keno kai allazei seira)
            }
            writer2.Flush();

        }
    }

}
