#!/bin/zsh

WORK_DIR="D:\TFG\Filtrado_noncoding"
OUTPUT_DIR="D:\TFG\Pytables_noncoding"

# Crear carpeta de salida si no existe
mkdir -p "$OUTPUT_DIR"

if [[ ! -d "$WORK_DIR" ]]; then
    echo "Error: El directorio '$WORK_DIR' no existe"
    exit 1
fi

cd "$WORK_DIR"

count=0

echo "Iniciando procesamiento de archivos PyTables..."
echo "Directorio origen: $(pwd)"
echo "Directorio destino: $OUTPUT_DIR"
echo "----------------------------------------"

for file in GSM*_out_filtered.h5; do
    if [[ ! -f "$file" ]]; then
        echo "No se encontraron archivos"
        continue
    fi

    GSM_NUM=""
    output_file=""

    # Extraer GSM
    if [[ "$file" =~ (GSM[0-9]+) ]]; then
        GSM_NUM=${BASH_REMATCH[1]}
    else
        echo "⚠️ No se pudo extraer GSM: $file"
        continue
    fi

    output_file="$OUTPUT_DIR/${GSM_NUM}.h5"

    echo "Procesando: $file"
    echo "  → $output_file"

    rm -f "$output_file"

    if ptrepack --complevel 5 --overwrite-nodes "${file}:/matrix" "${output_file}:/matrix"; then
        echo "OK"
        ((count++))
    else
        echo "Error en $file"
    fi

    echo "----------------------------------------"
done

echo "Procesamiento completado: $count archivos"
echo "Archivos guardados en: $OUTPUT_DIR"


