import xml.etree.ElementTree as ET
import os

def xml_to_fasta(xml_filename, output_filename=None):
    """
    Convierte un XML de resultados BLASTP a formato FASTA.
    
    Parámetros:
    - xml_filename: Nombre del archivo XML (ej: "resultados.xml")
    - output_filename: Nombre del archivo FASTA de salida (opcional, ej: "secuencias.fasta")
    """
    
    # Obtener la ruta del directorio donde está este script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Construir rutas completas para los archivos
    xml_path = os.path.join(script_dir, xml_filename)
    
    # Parsear el XML
    try:
        tree = ET.parse(xml_path)
        root = tree.getroot()
        print(f"Archivo XML leído correctamente: {xml_path}")
    except ET.ParseError as e:
        print(f"Error parsing XML: {e}")
        return
    except FileNotFoundError:
        print(f"Archivo no encontrado: {xml_path}")
        return
    
    # Generar nombre de salida si no se proporciona
    if output_filename is None:
        base_name = os.path.splitext(xml_filename)[0]
        output_filename = base_name + '.fasta'
    
    # Ruta completa para el archivo de salida
    output_path = os.path.join(script_dir, output_filename)
    
    sequences = []
    
    # Buscar todos los hits en el XML
    hits = root.findall(".//Hit")
    
    if not hits:
        print("No se encontraron hits en el XML")
        return
    
    print(f"Procesando {len(hits)} hits...")
    
    # Procesar cada hit
    for hit in hits:
        rank = hit.get("rank", "Unknown")
        title = hit.find("Title")
        accession = hit.find("Accession")
        species = hit.find("Species")
        
        # Obtener secuencias del alignment
        sequences_elem = hit.find("AlignmentSequences")
        if sequences_elem is None:
            continue
            
        subject_seq = sequences_elem.find("SubjectSequence")
        
        # Obtener estadísticas para el header
        stats = hit.find("Statistics")
        evalue = stats.find("EValue").text if stats and stats.find("EValue") is not None else "N/A"
        identity_pct = stats.find("IdentityPercentage").text if stats and stats.find("IdentityPercentage") is not None else "N/A"
        
        # Crear header informativo
        base_info = f"rank={rank} E-value={evalue} Identity={identity_pct}%"
        if species is not None and species.text:
            base_info += f" [{species.text}]"
        
        # Agregar secuencia subject (limpiar gaps)
        if subject_seq is not None and subject_seq.text:
            clean_subject = subject_seq.text.replace('-', '')
            
            if clean_subject:  # Solo agregar si la secuencia no está vacía
                subject_id = accession.text if accession is not None and accession.text != "N/A" else f"Hit_{rank}"
                subject_desc = title.text if title is not None else "No description"
                
                header = f">{subject_id} {subject_desc} | {base_info}"
                sequences.append((header, clean_subject))
    
    # Escribir archivo FASTA
    try:
        with open(output_path, 'w') as fasta_file:
            for header, sequence in sequences:
                fasta_file.write(f"{header}\n")
                # Escribir secuencia en líneas de 80 caracteres
                for i in range(0, len(sequence), 80):
                    fasta_file.write(f"{sequence[i:i+80]}\n")
        
        print(f"✓ Archivo FASTA creado exitosamente: {output_path}")
        print(f"✓ Total de secuencias extraídas: {len(sequences)}")
        
    except Exception as e:
        print(f"Error escribiendo archivo FASTA: {e}")


# CONFIGURACIÓN SIMPLE - SOLO CAMBIA ESTA LÍNEA:
xml_filename = "top10_albumina_hits_20250922_113448.xml"  # ← Pon aquí el nombre de tu archivo XML

# Ejecutar la conversión
if __name__ == "__main__":
    print("=== Conversor XML BLASTP a FASTA ===")
    print(f"Directorio de trabajo: {os.path.dirname(os.path.abspath(__file__))}")
    xml_to_fasta(xml_filename)
    print("=== Proceso completado ===")