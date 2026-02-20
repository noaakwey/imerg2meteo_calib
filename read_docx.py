import docx
import sys

def read_docx(file_path, out_file):
    try:
        doc = docx.Document(file_path)
        with open(out_file, 'w', encoding='utf-8') as f:
            for i, para in enumerate(doc.paragraphs):
                if para.text.strip():
                    f.write(f"[{i}] {para.text}\n")
        print("Success")
    except Exception as e:
        print(f"Error reading docx: {e}")

if __name__ == "__main__":
    if len(sys.argv) > 2:
        read_docx(sys.argv[1], sys.argv[2])
    else:
        print("Provide docx path and output path")
