import os
import argparse
from mistralai import Mistral

def process_pdfs_to_markdown(input_folder, output_folder, api_key):
    """
    Process all PDF files in the input folder, perform OCR on them,
    and save the results as markdown files in the output folder.
    """
    # Initialize the Mistral client
    client = Mistral(api_key=api_key)
    
    # Create output folder if it doesn't exist
    os.makedirs(output_folder, exist_ok=True)
    
    # Get all PDF files in the input folder
    pdf_files = [f for f in os.listdir(input_folder) if f.lower().endswith('.pdf')]
    
    if not pdf_files:
        print(f"No PDF files found in {input_folder}")
        return
    
    # Process each PDF file
    for pdf_file in pdf_files:
        input_path = os.path.join(input_folder, pdf_file)
        output_filename = os.path.splitext(pdf_file)[0] + '.md'
        output_path = os.path.join(output_folder, output_filename)
        
        print(f"Processing {pdf_file}...")
        
        try:
            # Upload the PDF file to Mistral API
            with open(input_path, "rb") as file_content:
                uploaded_pdf = client.files.upload(
                    file={
                        "file_name": pdf_file,
                        "content": file_content,
                    },
                    purpose="ocr"
                )
            
            # Generate a signed URL for the uploaded file
            signed_url = client.files.get_signed_url(file_id=uploaded_pdf.id)
            
            # Process the uploaded document for OCR
            ocr_response = client.ocr.process(
                model="mistral-ocr-latest",
                document={
                    "type": "document_url",
                    "document_url": signed_url.url,
                }
            )
            
            # Create a markdown document by combining content from all pages
            markdown_content = []
            
            # Add document title
            document_title = os.path.splitext(pdf_file)[0]
            markdown_content.append(f"# {document_title}\n\n")
            
            # Iterate over the pages and extract the markdown content
            for i, page in enumerate(ocr_response.pages, 1):
                # markdown_content.append(f"## Page {i}\n\n") # Uncomment to add page headers
                markdown_content.append(page.markdown + "\n\n")
            
            # Write to output file
            with open(output_path, 'w', encoding='utf-8') as md_file:
                md_file.write("".join(markdown_content))
            
            print(f"Successfully processed {pdf_file} -> {output_filename}")
            
        except Exception as e:
            print(f"Error processing {pdf_file}: {str(e)}")
    
    print(f"Processing complete. Results saved to {output_folder}")

def main():
    parser = argparse.ArgumentParser(description='Convert PDF files to Markdown using Mistral OCR')
    parser.add_argument('--input', required=True, help='Input folder containing PDF files')
    parser.add_argument('--output', required=True, help='Output folder for markdown files')
    parser.add_argument('--api-key', required=True, help='Mistral API key')
    
    args = parser.parse_args()
    
    process_pdfs_to_markdown(args.input, args.output, args.api_key)

if __name__ == "__main__":
    main()