import pandas as pd
import streamlit as st
# import utils from the same directory
# from utils.utils import check_variant_valid, randomize_features
import numpy as np
import pickle as pkl


def randomize_features():
    """This function will be deleted when the feature extraction function is implemented."""
    random_features_dictionary = {"blosum": np.random.randint(0, 100, 1)[0],
    "plddt_residue": np.random.randint(0, 100, 1)[0],
    "stability_delta": np.random.randint(0, 100, 1)[0],
    "hydrophobicity_delta": np.random.randint(0, 100, 1)[0],
    "volume_delta": np.random.randint(0, 100, 1)[0],
    "RSA_WT": np.random.randint(0, 100, 1)[0],
    "oda_delta": np.random.randint(0, 100, 1)[0],
    "sasa_delta": np.random.randint(0, 100, 1)[0],
    "pssm": np.random.randint(0, 100, 1)[0],
    "entropy": np.random.randint(0, 100, 1)[0],
    "secondary_structure_Beta": 1,
    "secondary_structure_alpha": 0,
    "secondary_structure_turn": 0,
    "secondary_structure_loop": 0}
    return random_features_dictionary


def check_variant_valid(variant, gene):
    """
    Check if the variant is valid for the given gene.
    Args:
        variant (str): The variant to check.
        gene (str): The gene to check.
    Returns:
        bool: True if the variant is valid, False otherwise
    """
    # Check if the variant is in the correct format
    if not variant.isalnum():
        return False
    # Check if the amino acids are valid
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"
    if variant[0] not in amino_acids or variant[-1] not in amino_acids:
        return False
    # Check if the position is valid (all characters except the first and last)
    if not variant[1:-1].isdigit():
        return False

    # delete later
    seq_myo7a = "MVILQQGDHVWMDLRLGQEFDVPIGAVVKLCDSGQVQVVDDEDNEHWISPQNATHIKPMHPTSVHGVEDMIRLGDLNEAGILRNLLIRYRDHLIYTYTGSILVAVNPYQLLSIYSPEHIRQYTNKKIGEMPPHIFAIADNCYFNMKRNSRDQCCIISGESGAGKTESTKLILQFLAAISGQHSWIEQQVLEATPILEAFGNAKTIRNDNSSRFGKYIDIHFNKRGAIEGAKIEQYLLEKSRVCRQALDERNYHVFYCMLEGMSEDQKKKLGLGQASDYNYLAMGNCITCEGRVDSQEYANIRSAMKVLMFTDTENWEISKLLAAILHLGNLQYEARTFENLDACEVLFSPSLATAASLLEVNPPDLMSCLTSRTLITRGETVSTPLSREQALDVRDAFVKGIYGRLFVWIVDKINAAIYKPPSQDVKNSRRSIGLLDIFGFENFAVNSFEQLCINFANEHLQQFFVRHVFKLEQEEYDLESIDWLHIEFTDNQDALDMIANKPMNIISLIDEESKFPKGTDTTMLHKLNSQHKLNANYIPPKNNHETQFGINHFAGIVYYETQGFLEKNRDTLHGDIIQLVHSSRNKFIKQIFQADVAMGAETRKRSPTLSSQFKRSLELLMRTLGACQPFFVRCIKPNEFKKPMLFDRHLCVRQLRYSGMMETIRIRRAGYPIRYSFVEFVERYRVLLPGVKPAYKQGDLRGTCQRMAEAVLGTHDDWQIGKTKIFLKDHHDMLLEVERDKAITDRVILLQKVIRGFKDRSNFLKLKNAATLIQRHWRGHNCRKNYGLMRLGFLRLQALHRSRKLHQQYRLARQRIIQFQARCRAYLVRKAFRHRLWAVLTVQAYARGMIARRLHQRLRAEYLWRLEAEKMRLAEEEKLRKEMSAKKAKEEAERKHQERLAQLAREDAERELKEKEAARRKKELLEQMERARHEPVNHSDMVDKMFGFLGTSGGLPGQEGQAPSGFEDLERGRREMVEEDLDAALPLPDEDEEDLSEYKFAKFAATYFQGTTTHSYTRRPLKQPLLYHDDEGDQLAALAVWITILRFMGDLPEPKYHTAMSDGSEKIPVMTKIYETLGKKTYKRELQALQGEGEAQLPEGQKKSSVRHKLVHLTLKKKSKLTEEVTKRLHDGESTVQGNSMLEDRPTSNLEKLHFIIGNGILRPALRDEIYCQISKQLTHNPSKSSYARGWILVSLCVGCFAPSEKFVKYLRNFIHGGPPGYAPYCEERLRRTFVNGTRTQPPSWLELQATKSKKPIMLPVTFMDGTTKTLLTDSATTAKELCNALADKISLKDRFGFSLYIALFDKVSSLGSGSDHVMDAISQCEQYAKEQGAQERNAPWRLFFRKEVFTPWHSPSEDNVATNLIYQQVVRGVKFGEYRCEKEDDLAELASQQYFVDYGSEMILERLLNLVPTYIPDREITPLKTLEKWAQLAIAAHKKGIYAQRRTDAQKVKEDVVSYARFKWPLLFSRFYEAYKFSGPSLPKNDVIVAVNWTGVYFVDEQEQVLLELSFPEIMAVSSSRECRVWLSLGCSDLGCAAPHSGWAGLTPAGPCSPCWSCRGAKTTAPSFTLATIKGDEYTFTSSNAEDIRDLVVTFLEGLRKRSKYVVALQDNPNPAGEESGFLSFAKGDLIILDHDTGEQVMNSGWANGINERTKQRGDFPTDSVYVMPTVTMPPREIVALVTMTPDQRQDVVRLLQLRTAEPEVRAKPYTLEEFSYDYFRPPPKHTLSRVMVSKARGKDRLWSHTREPLKQALLKKLLGSEELSQEACLAFIAVLKYMGDYPSKRTRSVNELTDQIFEGPLKAEPLKDEAYVQILKQLTDNHIRYSEERGWELLWLCTGLFPPSNILLPHVQRFLQSRKHCPLAIDCLQRLQKALRNGSRKYPPHLVEVEAIQHKTTQIFHKVYFPDDTDEAFEVESSTKAKDFCQNIATRLLLKSSEGFSLFVKIADKVLSVPENDFFFDFVRHLTDWIKKARPIKDGIVPSLTYQVFFMKKLWTTTVPGKDPMADSIFHYYQELPKYLRGYHKCTREEVLQLGALIYRVKFEEDKSYFPSIPKLLRELVPQDLIRQVSPDDWKRSIVAYFNKHAGKSKEEAKLAFLKLIFKWPTFGSAFFEVKQTTEPNFPEILLIAINKYGVSLIDPKTKDILTTHPFTKISNWSSGNTYFHITIGNLVRGSKLLCETSLGYKMDDLLTSYISQMLTAMSKQRGSRSGK"

    # Check if this position is valid for the given gene, using uniprot
    gene_sequence = seq_myo7a
    # Check if the position is within the length of the gene sequence
    if int(variant[1:-1]) > len(gene_sequence):
        return False
    # Check if the amino acid at this position is the same as the one in the variant
    if gene_sequence[int(variant[1:-1]) - 1] != variant[0]:
        return False

    return True


genes_list = ["SLC26A4", "GJB2", "COL2A1", "COL4A5", "COL4A3", "MYO7A", "WFS1", "FGFR1"]
features = ["blosum", "plddt_residue", "stability_delta", "hydrophobicity_delta", "volume_delta", "RSA_WT",
            "oda_delta", "sasa_delta", "pssm", "entropy", "secondary_structure_Beta", "secondary_structure_alpha",
            "secondary_structure_turn", "secondary_structure_loop"]

# Randomize features and store in a dictionary
random_features_dict = randomize_features()
# TODO: Later, include the feature extraction and preprocessing functions here


def predict():
    model_path = f"C:\\Users\\InbarBlech\\PycharmProjects\\Thesis\\predictions_vs_real\\gene_specific_predictors\\{gene}_model.pkl"
    model = pkl.load(open(model_path, "rb"))
    X = pd.DataFrame(random_features_dict, columns=features)
    predictions = model.predict(X)[0]
    proba = model.predict_proba(X)[0]

    if predictions == 1:
        st.progress(proba[1])
        st.write("The variant is predicted to be pathogenic.")
    else:
        st.progress(proba[0])
        st.write("The variant is predicted to be benign.")


with st.container():
    st.write("## Predictor")
    left_column, right_column = st.columns(2)
    with left_column:
        st.subheader("Enter a variant")
        variant = st.text_input("Variant", placeholder="e.g. R75Q")
    with right_column:
        st.subheader("Select a gene")
        gene = st.selectbox("Select a gene", genes_list)
    # Add a button
    if st.button("Predict"):
        if variant == "" or gene == "":
            st.write("Please enter a variant and a gene.")
        # Check that variant is in the correct format
        # if not check_variant_valid(variant, gene):
        #     st.write("Error. Please enter a valid variant.")
        if variant and gene:
            # Perform feature extraction
            st.write("Computing features for the variant...")
            # features = extract_features(variant, gene)  # TODO: Implement this function
            # Perform feature preprocessing
            st.write("Preprocessing the features...")
            # features = preprocess_features(features)  # TODO: Implement this function
            # Meanwhile, use random features
            st.write("Using random features for demonstration purposes.")
            # Make predictions
            st.write(f"Predicting the pathogenicity of the variant {variant} in gene {gene}...")
            predict()
        else:
            st.write("Please enter valid variant and gene.")
    st.write("---")
