using System.Collections.Generic;
using ISAAR.MSolve.FEM.Interfaces;
using ISAAR.MSolve.FEM.Entities;
using System.Linq;

namespace ISAAR.MSolve.SamplesConsole.SupportiveClasses
{
    public static class DdmCalculationsPartb_v2 // exei allaxei mono to onoma den exei ginei update se model_v2 xrhsh klp.
    {
        public static (Dictionary<int, Dictionary<int, IList<int>>>, Dictionary<int, List<int>>, Dictionary<int, List<int>>) FindEmbeddedElementsSubdomains(Model_v2 model,int totalSubdomains)
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

            foreach (Element_v2 element in model.ElementsDictionary.Values) // ean xeroume apo thn arxh to id ton embedded mporoume na ton dinoume
            {
                if (element.ElementType is IEmbeddedElement_v2)
                {
                    var e1 = element.ElementType as IEmbeddedElement_v2;
                    Dictionary<int, IList<int>> HostSubdomains = new Dictionary<int, IList<int>>();
                    foreach (var embeddedNode in (e1).EmbeddedNodes)
                    {
                        //1
                        Element_v2 hostELement = embeddedNode.EmbeddedInElement;
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
                        if ( AssignedSubdomains.ContainsKey(HostSubdomains.ElementAt(0).Key) )
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
            return (AmbiguousEmbeddedElementsHostSubdomainsAndSpecifcHexaElementsInThem, hexaConnectsShells,AssignedSubdomains);
        }

        public static (Dictionary<int, Dictionary<int, IList<int>>>, Dictionary<int, List<int>>, Dictionary<int, List<int>>) FindEmbeddedElementsSubdomainsCorrected(Model_v2 model, int totalSubdomains)
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

            foreach (Element_v2 element in model.ElementsDictionary.Values) // ean xeroume apo thn arxh to id ton embedded mporoume na ton dinoume
            {
                if (element.ElementType is IEmbeddedElement_v2)
                {
                    Dictionary<int, List<int>> hexaConnectsShellsLocal = new Dictionary<int, List<int>>();

                    var e1 = element.ElementType as IEmbeddedElement_v2;
                    Dictionary<int, IList<int>> HostSubdomains = new Dictionary<int, IList<int>>();
                    foreach (var embeddedNode in (e1).EmbeddedNodes)
                    {
                        //1
                        Element_v2 hostELement = embeddedNode.EmbeddedInElement;
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
                        foreach(int hexaID in hexaConnectsShellsLocal.Keys)
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
                foreach(int subdID in AmbigElements[ambElementID].Keys)
                {
                    //numSpecificElements[subdCounter] = AmbigElements[ambElementID][subdID].Count; // Count dld ta element ths inner listas tou AmbigElements...
                    numSpecificElements[subdCounter] = AmbigElements[ambElementID].ElementAt(subdCounter).Value.Count; // Value einai h inner lista tou AmbigElements...
                    subdCounter++;
                }

                int thesiChosenSubd = 0;
                int nElements = 0;
                for (int i1=0; i1<numSpecificElements.GetLength(0); i1++ )
                {
                    if (numSpecificElements[i1]>nElements)
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
            int[][] subdIdsAndElements = new int[maxSubdId+1][]; //todo:revisit this

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

        public static  Dictionary<int, List<int>> FindEmbeddedElementsSubdomainsCorrectedSimple(Model_v2 model, int totalSubdomains)
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

            foreach (Element_v2 element in model.ElementsDictionary.Values) // ean xeroume apo thn arxh to id ton embedded mporoume na ton dinoume
            {
                if (element.ElementType is IEmbeddedElement_v2)
                {
                    Dictionary<int, List<int>> hexaConnectsShellsLocal = new Dictionary<int, List<int>>();

                    var e1 = element.ElementType as IEmbeddedElement_v2;
                    Dictionary<int, IList<int>> HostSubdomains = new Dictionary<int, IList<int>>();
                    foreach (var embeddedNode in (e1).EmbeddedNodes)
                    {
                        //1
                        Element_v2 hostELement = embeddedNode.EmbeddedInElement;
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
                            if(HostSubdomains[subdId].Count>hexaListlength)
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
}
