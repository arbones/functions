### `pdf2md.py` : PDF to Markdown Converter using Mistral OCR API

This function processes a folder of PDF files <img src="https://upload.wikimedia.org/wikipedia/commons/thumb/8/87/PDF_file_icon.svg/195px-PDF_file_icon.svg.png" width="15">, performs Optical Character Recognition (OCR) on each document using the Mistral AI OCR API <img src="https://upload.wikimedia.org/wikipedia/commons/thumb/e/e6/Mistral_AI_logo_%282025%E2%80%93%29.svg/1024px-Mistral_AI_logo_%282025%E2%80%93%29.svg.png" width="30">, and saves the extracted text as markdown files <img src="https://upload.wikimedia.org/wikipedia/commons/thumb/4/41/1280px_Markdown_with_White_Background.png/320px-1280px_Markdown_with_White_Background.png" width="40"> in an output folder.

#### Arguments

- `input_folder` (string): Directory path containing PDF files to be processed
- `output_folder` (string): Directory path where the resulting markdown files will be saved
- `api_key` (string): Your Mistral AI API key for authentication

#### Description

The function:
1. Initializes the Mistral API client with your API key
2. Creates the output folder if it doesn't exist
3. Processes each PDF file in the input folder:
   - Uploads the PDF to Mistral's servers
   - Generates a signed URL for OCR processing
   - Extracts text from each page using OCR
   - Combines the text into a well-formatted markdown document
   - Saves the result to the output folder with the same filename but .md extension
4. Provides progress feedback during processing

#### Example Usage

```bash
python pdf_to_markdown.py --input /path/to/pdf/folder --output /path/to/markdown/folder --api-key your_mistral_api_key
```

#### Output Format

For each PDF file processed, the script creates a markdown file with:
- Document title (derived from the PDF filename)
- Numbered sections for each page
- Text content from each page formatted in markdown

#### Requirements

- Python 3.6+
- Mistral AI Python client: `pip install mistralai`
- Valid Mistral AI API key with OCR capabilities

#### Error Handling

The script includes error handling to continue processing other files if one file fails. Errors are reported to the console but don't stop the entire batch from being processed.