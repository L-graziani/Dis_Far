from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO
from io import StringIO
import time
import xml.etree.ElementTree as ET
from xml.dom import minidom
from pathlib import Path
# Aquí va la secuencia de la albúmina humana en formato FASTA
albumina_fasta = """
>sp|P02768|ALBU_HUMAN Albumin OS=Homo sapiens OX=9606 GN=ALB PE=1 SV=2
MKWVTFISLLFLFSSAYSRGVFRRDAHKSEVAHRFKDLGEENFKALVLIAFAQYLQQCPF
EDHVKLVNEVTEFAKTCVADESAENCDKSLHTLFGDKLCTVATLRETYGEMADCCAKQEP
ERNECFLQHKDDNPNLPRLVRPEVDVMCTAFHDNEETFLKKYLYEIARRHPYFYAPELLF
FAKRYKAAFTECCQAADKAACLLPKLDELRDEGKASSAKQRLKCASLQKFGERAFKAWAV
ARLSQRFPKAEFAEVSKLVTDLTKVHTECCHGDLLECADDRADLAKYICENQDSISSKLK
ECCEKPLLEKSHCIAEVENDEMPADLPSLAADFVESKDVCKNYAEAKDVFLGMFLYEYAR
RHPDYSVVLLLRLAKTYETTLEKCCAAADPHECYAKVFDEFKPLVEEPQNLIKQNCELFE
QLGEYKFQNALLVRYTKKVPQVSTPTLVEVSRNLGKVGSKCCKHPEAKRMPCAEDYLSVV
LNQLCVLHEKTPVSDRVTKCCTESLVNRRPCFSALEVDETYVPKEFNAETFTFHADICTL
SEKERQIKKQTALVELVKHKPKATKEQLKAVMDDFAAFVEKCCKADDKETCFAEEGKKLV
AASQAALGL
"""

def guardar_top10_xml(top_alignments, query_id):
    """
    Guarda las 10 mejores secuencias en formato XML estructurado
    """
    # Crear elemento raíz
    root = ET.Element("BlastResults")
    root.set("query_id", query_id)
    root.set("timestamp", time.strftime("%Y-%m-%d %H:%M:%S"))
    
    # Elemento para metadata
    metadata = ET.SubElement(root, "Metadata")
    ET.SubElement(metadata, "Program").text = "BLASTP"
    ET.SubElement(metadata, "Database").text = "RefSeq Proteins (Vertebrates)"
    ET.SubElement(metadata, "TopHits").text = str(len(top_alignments))
    
    # Elemento para los hits
    hits_element = ET.SubElement(root, "TopHits")
    
    for i, alignment in enumerate(top_alignments):
        hsp = alignment.hsps[0]  # Mejor HSP
        
        # Crear elemento hit
        hit = ET.SubElement(hits_element, "Hit")
        hit.set("rank", str(i + 1))
        
        # Información básica del hit
        ET.SubElement(hit, "Title").text = alignment.title
        ET.SubElement(hit, "Length").text = str(alignment.length)
        ET.SubElement(hit, "Accession").text = alignment.accession if hasattr(alignment, 'accession') else "N/A"
        
        # Estadísticas del HSP
        stats = ET.SubElement(hit, "Statistics")
        ET.SubElement(stats, "Score").text = str(hsp.score)
        ET.SubElement(stats, "BitScore").text = str(hsp.bits)
        ET.SubElement(stats, "EValue").text = str(hsp.expect)
        ET.SubElement(stats, "Identities").text = f"{hsp.identities}/{hsp.align_length}"
        ET.SubElement(stats, "IdentityPercentage").text = f"{hsp.identities/hsp.align_length*100:.2f}"
        ET.SubElement(stats, "Positives").text = f"{hsp.positives}/{hsp.align_length}"
        ET.SubElement(stats, "PositivePercentage").text = f"{hsp.positives/hsp.align_length*100:.2f}"
        ET.SubElement(stats, "Gaps").text = f"{hsp.gaps}/{hsp.align_length}"
        ET.SubElement(stats, "GapPercentage").text = f"{hsp.gaps/hsp.align_length*100:.2f}"
        
        # Posiciones del alignment
        positions = ET.SubElement(hit, "AlignmentPositions")
        ET.SubElement(positions, "QueryStart").text = str(hsp.query_start)
        ET.SubElement(positions, "QueryEnd").text = str(hsp.query_end)
        ET.SubElement(positions, "SubjectStart").text = str(hsp.sbjct_start)
        ET.SubElement(positions, "SubjectEnd").text = str(hsp.sbjct_end)
        
        # Secuencias del alignment
        sequences = ET.SubElement(hit, "AlignmentSequences")
        ET.SubElement(sequences, "QuerySequence").text = hsp.query
        ET.SubElement(sequences, "MatchString").text = hsp.match
        ET.SubElement(sequences, "SubjectSequence").text = hsp.sbjct
        
        # Extraer información de la especie si está disponible
        if '[' in alignment.title and ']' in alignment.title:
            species = alignment.title[alignment.title.rfind('[')+1:alignment.title.rfind(']')]
            ET.SubElement(hit, "Species").text = species
    
    # Formatear y guardar XML
    xml_string = ET.tostring(root, encoding='unicode')
    dom = minidom.parseString(xml_string)
    pretty_xml = dom.toprettyxml(indent="  ")
    
    # Guardar en el mismo directorio del script (no en la raíz del workspace)
    output_dir = Path(__file__).resolve().parent if '__file__' in globals() else Path.cwd()
    filename = output_dir / f"top10_albumina_hits_{time.strftime('%Y%m%d_%H%M%S')}.xml"

    with open(filename, 'w', encoding='utf-8') as f:
        f.write(pretty_xml)

    print(f"\n✓ Top 10 secuencias guardadas en: {filename}")
    return str(filename)

def validar_y_limpiar_fasta(secuencia_fasta):
    """
    Valida y limpia una secuencia FASTA
    """
    # Limpiar líneas vacías y espacios
    lineas_limpias = [line.strip() for line in secuencia_fasta.strip().split('\n') if line.strip()]
    
    if not lineas_limpias:
        raise ValueError("La secuencia FASTA está vacía")
    
    # Verificar que empiece con '>'
    if not lineas_limpias[0].startswith('>'):
        raise ValueError("La primera línea debe empezar con '>' (header FASTA)")
    
    # Verificar que tenga al menos una línea de secuencia
    if len(lineas_limpias) < 2:
        raise ValueError("Falta la secuencia después del header")
    
    return '\n'.join(lineas_limpias)

def realizar_blastp_albumina():
    """
    Realiza un BLAST-P de la albúmina humana contra la base de datos de proteínas de vertebrados
    """
    
    # Parsear la secuencia FASTA
    try:
        # Validar y limpiar la secuencia
        fasta_limpia = validar_y_limpiar_fasta(albumina_fasta)
        
        seq_record = SeqIO.read(StringIO(fasta_limpia), "fasta")
        print(f"Secuencia cargada: {seq_record.id}")
        print(f"Descripción: {seq_record.description}")
        print(f"Longitud: {len(seq_record.seq)} aminoácidos")
        print("="*50)
    except Exception as e:
        print(f"Error al parsear la secuencia: {e}")
        print("\nFormato esperado:")
        print(">sp|P02768|ALBU_HUMAN Serum albumin - Homo sapiens (Human)")
        print("MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIAFSQYLQQCPFDEHVKLVNEL...")
        print("\nVerifica que tu variable 'albumina_fasta' tenga este formato.")
        return
    
    # Configurar parámetros del BLAST
    print("Iniciando BLAST-P...")
    print("Base de datos: refseq_protein (vertebrados)")
    print("Esto puede tomar varios minutos...")
    
    try:
        # Realizar BLAST-P
        result_handle = NCBIWWW.qblast(
            program="blastp",           # Programa BLAST para proteínas
            database="refseq_protein",  # Base de datos RefSeq de proteínas
            sequence=str(seq_record.seq),  # Secuencia query
            entrez_query="vertebrates[organism]",  # Filtrar solo vertebrados
            hitlist_size=50,            # Número máximo de hits
            expect=0.001,               # E-value threshold
            word_size=3,                # Tamaño de palabra para proteínas
            matrix_name="BLOSUM62"      # Matriz de sustitución
        )
        
        # Parsear resultados
        blast_records = NCBIXML.parse(result_handle)
        blast_record = next(blast_records)
        
        print(f"\nResultados del BLAST-P:")
        print(f"Query: {blast_record.query}")
        print(f"Longitud query: {blast_record.query_length}")
        print(f"Base de datos: {blast_record.database}")
        print(f"Número de secuencias en DB: {blast_record.database_sequences}")
        print("="*80)
        
        # Procesar alignments
        if blast_record.alignments:
            print(f"Se encontraron {len(blast_record.alignments)} alignments")
            print("\nTop 10 resultados:")
            print("-"*80)
            
            # Guardar las 10 mejores secuencias en XML
            guardar_top10_xml(blast_record.alignments[:10], seq_record.id)
            
            for i, alignment in enumerate(blast_record.alignments[:10]):
                # Tomar el primer HSP (High-scoring Segment Pair)
                hsp = alignment.hsps[0]
                
                print(f"\n{i+1}. {alignment.title[:80]}...")
                print(f"   Longitud: {alignment.length} aa")
                print(f"   E-value: {hsp.expect:.2e}")
                print(f"   Score: {hsp.score}")
                print(f"   Identidades: {hsp.identities}/{hsp.align_length} ({hsp.identities/hsp.align_length*100:.1f}%)")
                print(f"   Positivos: {hsp.positives}/{hsp.align_length} ({hsp.positives/hsp.align_length*100:.1f}%)")
                print(f"   Gaps: {hsp.gaps}/{hsp.align_length} ({hsp.gaps/hsp.align_length*100:.1f}%)")
                
                # Mostrar una parte del alignment
                if len(hsp.query) > 60:
                    print(f"   Query:  {hsp.query[:60]}...")
                    print(f"   Match:  {hsp.match[:60]}...")
                    print(f"   Sbjct:  {hsp.sbjct[:60]}...")
                else:
                    print(f"   Query:  {hsp.query}")
                    print(f"   Match:  {hsp.match}")
                    print(f"   Sbjct:  {hsp.sbjct}")
        
        else:
            print("No se encontraron alignments significativos")
            
        result_handle.close()
        
    except Exception as e:
        print(f"Error durante el BLAST: {e}")
        return


# Ejecutar el análisis principal
if __name__ == "__main__":
    print("BLAST-P: Albúmina humana vs proteínas de vertebrados")
    print("="*60)
    
    # Opción 1: Análisis directo
    realizar_blastp_albumina()
    
