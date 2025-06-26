import gleam/int
import gleam/fetch
import gleam/http/request
import gleam/http/response
import gleam/io
import gleam/javascript/promise.{type Promise}
import gleam/list
import gleam/result
import gsv

type Complex {
  Complex(id: String, name: String)
}

type Error {
  CSVError(gsv.ParseError)
  ParsingError(String)
  FetchError(fetch.FetchError)
}

fn strerror(error: Error) -> String {
  case error {
    CSVError(csv_error) -> {
      case csv_error {
        gsv.UnescapedQuote(position) -> "CSV: Unescaped quote @ " <> int.to_string(position)
        gsv.UnclosedEscapedField(start) -> "CSV: Unclosed escape starting " <> int.to_string(start)
      }
    }
    ParsingError(error_text) -> "Parsing Error: " <> error_text
    FetchError(fetch_error) -> {
      case fetch_error {
        fetch.NetworkError(e) -> "Error: " <> e
        fetch.UnableToReadBody -> "Error: unable to read body"
        fetch.InvalidJsonBody
      }
    }
  }
}

fn map_to_error(
  prom: Promise(Result(a, fetch.FetchError)),
) -> Promise(Result(a, Error)) {
  promise.map(prom, fn(res) { result.map_error(res, FetchError) })
}

fn fetch_complexome(taxon_id: String) -> Promise(Result(List(Complex), Error)) {
  let assert Ok(request) =
    request.to(
      "https://ftp.ebi.ac.uk/pub/databases/intact/complex/current/complextab/"
      <> taxon_id
      <> ".tsv",
    )
  fetch.send(request)
  |> map_to_error
  |> promise.try_await(fn(resp) { map_to_error(fetch.read_text_body(resp)) })
  |> promise.map(fn(body) {
    case body {
      Ok(response.Response(200, _, text)) -> parse_complexome(text)
      _ -> Error(ParsingError("What even is this?"))
    }
  })
}

fn parse_complexome(data: String) -> Result(List(Complex), Error) {
  use rows <- result.try(gsv.to_lists(data) |> result.map_error(CSVError))
  list.try_map(rows, fn(rows) {
    case rows {
      [id, name, _, taxon_id, participants, _, go_terms, ext_participants] ->
        Ok(Complex(id, name))
      _ -> Error(ParsingError("Error"))
    }
  })
}

pub fn main() -> Nil {
  promise.tap(fetch_complexome("9606"), fn(res) {
    case res {
      Ok(complexes) -> {
        echo complexes
        Nil
      }
      Error(error) -> io.println_error(strerror(error))
    }
  })
  Nil
}
