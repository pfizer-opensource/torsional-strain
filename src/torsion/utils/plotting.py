from openeye import oechem, oedepict
from .process_sd_data import get_sd_data
from ..core import get_dihedral

def plot_indices(mol2, width=200, height=200):
    mol = mol2.CreateCopy()

    opts = oedepict.OE2DMolDisplayOptions(width, height, oedepict.OEScale_AutoScale)
    opts.SetAtomPropertyFunctor(oedepict.OEDisplayAtomIdx())
    oedepict.OEPrepareDepiction(mol)

    disp = oedepict.OE2DMolDisplay(mol, opts)
    img = oedepict.OEImage(width, height)
    
    oedepict.OERenderMolecule(img, disp)
    return(img)

def plot_dihedral(mol2, width=200, height=200):
    mol = mol2.CreateCopy()
    dihedralAtomIndices = [int(x)-1 for x in get_sd_data(mol, 'TORSION_ATOMS_FRAGMENT').split()]
    dih, tor = get_dihedral(mol, dihedralAtomIndices)

    opts = oedepict.OE2DMolDisplayOptions(width, height, oedepict.OEScale_AutoScale)
    oedepict.OEPrepareDepiction(mol)

    disp = oedepict.OE2DMolDisplay(mol, opts)
    img = oedepict.OEImage(width, height)

    hstyle = oedepict.OEHighlightByBallAndStick(oechem.OEBlueTint)
    oedepict.OEAddHighlighting(disp, hstyle, dih)
    hstyle = oedepict.OEHighlightByColor(oechem.OERed)
    oedepict.OEAddHighlighting(disp, hstyle, tor)

    oedepict.OERenderMolecule(img, disp)
    return(img)

def highlight_atoms_in_mol(mol2, dihedralAtomIndices, width=200, height=200):
    mol = mol2.CreateCopy()
    opts = oedepict.OE2DMolDisplayOptions(width, height, oedepict.OEScale_AutoScale)
    oedepict.OEPrepareDepiction(mol)

    disp = oedepict.OE2DMolDisplay(mol, opts)
    img = oedepict.OEImage(width, height)

    hstyle = oedepict.OEHighlightByBallAndStick(oechem.OEBlueTint)
    for atom_idx in dihedralAtomIndices:
        oedepict.OEAddHighlighting(disp, hstyle, oechem.OEHasAtomIdx(atom_idx))

    oedepict.OERenderMolecule(img, disp)
    return(img)
   
    
def draw_subsearch_highlights(mol, subsearch, width=400., height=400.):
    """
    Draws the hits for the substructure in a given molecule.
    
    Copied from http://notebooks.eyesopen.com/substructure-search-pandas-oenotebook.html
    """
    opts = oedepict.OE2DMolDisplayOptions(width, height, oedepict.OEScale_AutoScale)

    mol = oechem.OEGraphMol(mol)
    oedepict.OEPrepareDepiction(mol)
    img = oedepict.OEImage(width, height)
    hstyle = oedepict.OEHighlightByBallAndStick(oechem.OEBlueTint)

    disp = oedepict.OE2DMolDisplay(mol, opts)
    unique = True
    for match in subsearch.Match(mol, unique):
        oedepict.OEAddHighlighting(disp, hstyle, match)
    
    oedepict.OERenderMolecule(img, disp)
    #return oenb.draw_oeimage_to_img_tag(img)
    return img
